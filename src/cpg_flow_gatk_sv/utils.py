"""
Common methods for all GATK-SV workflows
"""

import functools
import itertools
import re
from collections import defaultdict
from enum import Enum
from functools import cache
from os.path import join
from random import randint
from typing import TYPE_CHECKING, Any

import loguru

from cpg_flow import targets, utils, workflow
from cpg_utils import Path, config, cromwell, hail_batch, to_path
from metamist.graphql import gql, query

if TYPE_CHECKING:
    from hailtop.batch.job import BashJob

GATK_SV_COMMIT = config.config_retrieve(['workflow', 'gatk_sv_commit'])
SV_CALLERS = ['manta', 'wham', 'scramble']

_FASTA_STRING = None
PED_FAMILY_ID_REGEX = re.compile(r'(^[A-Za-z0-9_]+$)')


VCF_QUERY = gql(
    """
    query MyQuery($dataset: String!) {
        project(name: $dataset) {
            analyses(active: {eq: true}, type: {eq: "sv"}, status: {eq: COMPLETED}) {
                output
                timestampCompleted
            }
        }
    }
""",
)


class CromwellJobSizes(Enum):
    """
    Enum for polling intervals
    """

    SMALL = 'small'
    MEDIUM = 'medium'
    LARGE = 'large'


@cache
def get_sv_callers(add_manta: bool = False):
    if only_jobs := config.config_retrieve(['workflow', 'GatherSampleEvidence', 'only_jobs'], None):
        callers = [caller for caller in SV_CALLERS if caller in only_jobs]
        if 'scramble' in only_jobs and 'manta' not in only_jobs and add_manta:
            callers.append('manta')
        if not callers:
            loguru.logger.warning('No SV callers enabled')
        return callers
    return SV_CALLERS


@cache
def create_polling_intervals() -> dict:
    """
    Set polling intervals for cromwell status
    these values are integers, indicating seconds
    for each job size, there is a min and max value
    analysis-runner implements a backoff-retrier when checking for
    success, with a minimum value, gradually reaching a max ceiling

    a config section containing overrides would look like

    [cromwell_polling_intervals.medium]
    min = 69
    max = 420
    """

    # create this dict with default values
    polling_interval_dict = {
        CromwellJobSizes.SMALL: {'min': 30, 'max': 140},
        CromwellJobSizes.MEDIUM: {'min': 40, 'max': 400},
        CromwellJobSizes.LARGE: {'min': 200, 'max': 2000},
    }

    # update if these exist in config
    for job_size in CromwellJobSizes:
        if val := config.config_retrieve(['cromwell_polling_intervals', job_size.value], False):
            polling_interval_dict[job_size].update(val)
    return polling_interval_dict


def get_fasta_string() -> Path:
    """
    find or return the fasta to use
    """
    global _FASTA_STRING
    if _FASTA_STRING is None:
        _FASTA_STRING = config.config_retrieve(['workflow', 'ref_fasta'])
    return _FASTA_STRING


def get_references(keys: list[str | dict[str, str]]) -> dict[str, str | list[str]]:
    """
    Dict of WDL inputs with reference file paths.
    """
    res: dict[str, str | list[str]] = {}
    for key in keys:
        # Keys can be maps (e.g. {'MakeCohortVcf.cytobands': 'cytoband'})
        use_key, ref_key = next(iter(key.items())) if isinstance(key, dict) else (key, key)

        # e.g. GATKSVPipelineBatch.rmsk -> rmsk
        ref_key = ref_key.split('.')[-1]

        try:
            res[use_key] = config.reference_path(f'gatk_sv/{ref_key}')  # type: ignore[index]
        except (KeyError, config.ConfigError):
            res[use_key] = config.reference_path(f'broad/{ref_key}')  # type: ignore[index]

    return res


def add_gatk_sv_jobs(
    dataset: targets.Dataset,
    wfl_name: str,
    # "dict" is invariant (supports updating), "Mapping" is covariant (read-only)
    # we have to support inputs of type dict[str, str], so using Mapping here:
    input_dict: dict[str, Any],
    expected_out_dict: dict[str, Path | list[Path]],
    sequencing_group_id: str | None = None,
    labels: dict[str, str] | None = None,
    job_size: CromwellJobSizes = CromwellJobSizes.MEDIUM,
) -> list['BashJob']:
    """
    Generic function to add a job that would run one GATK-SV workflow.
    """

    # create/retrieve dictionary of polling intervals for cromwell status
    polling_intervals = create_polling_intervals()

    # obtain upper and lower polling bounds for this job size
    polling_minimum = randint(polling_intervals[job_size]['min'], polling_intervals[job_size]['min'] * 2)
    polling_maximum = randint(polling_intervals[job_size]['max'], polling_intervals[job_size]['max'] * 2)

    # If a config section exists for this workflow, apply overrides
    if override := config.config_retrieve(['resource_overrides', wfl_name], False):
        input_dict |= override

    # Where Cromwell writes the output.
    # Will be different from paths in expected_out_dict:
    output_prefix = f'gatk_sv/output/{wfl_name}/{dataset.name}'
    if sequencing_group_id:
        output_prefix = join(output_prefix, sequencing_group_id)

    outputs_to_collect = {}
    for key, value in expected_out_dict.items():
        if isinstance(value, list):
            outputs_to_collect[key] = cromwell.CromwellOutputType.array_path(
                name=f'{wfl_name}.{key}', length=len(value)
            )
        else:
            outputs_to_collect[key] = cromwell.CromwellOutputType.single_path(f'{wfl_name}.{key}')

    # pre-process input_dict
    paths_as_strings: dict = {}
    for key, value in input_dict.items():
        if isinstance(value, Path):
            paths_as_strings[f'{wfl_name}.{key}'] = str(value)
        elif isinstance(value, list | set):
            paths_as_strings[f'{wfl_name}.{key}'] = [str(v) for v in value]
        else:
            paths_as_strings[f'{wfl_name}.{key}'] = value

    job_prefix = utils.make_job_name(wfl_name, sequencing_group=sequencing_group_id, dataset=dataset.name)

    submit_j, output_dict = cromwell.run_cromwell_workflow_from_repo_and_get_outputs(
        b=hail_batch.get_batch(),
        job_prefix=job_prefix,
        dataset=config.config_retrieve(['workflow', 'dataset']),
        repo='gatk-sv',
        commit=GATK_SV_COMMIT,
        cwd='wdl',
        workflow=f'{wfl_name}.wdl',
        libs=['.'],
        output_prefix=output_prefix,
        input_dict=paths_as_strings,
        outputs_to_collect=outputs_to_collect,
        driver_image=config.config_retrieve(['workflow', 'driver_image']),
        copy_outputs_to_gcp=config.config_retrieve(['workflow', 'copy_outputs'], False),
        labels=labels,
        min_watch_poll_interval=polling_minimum,
        max_watch_poll_interval=polling_maximum,
        time_limit_seconds=config.config_retrieve(['workflow', 'time_limit_seconds'], None),
    )

    copy_j = hail_batch.get_batch().new_bash_job(f'{job_prefix}: copy outputs')
    copy_j.image(config.config_retrieve(['workflow', 'driver_image']))
    cmds = []
    for key, resource in output_dict.items():
        out_path = expected_out_dict[key]
        if isinstance(resource, list):
            for source, dest in zip(resource, out_path, strict=False):
                cmds.append(f'gcloud storage cp "$(cat {source})" "{dest}"')
        else:
            cmds.append(f'gcloud storage cp "$(cat {resource})" "{out_path}"')
    copy_j.command(hail_batch.command(cmds, setup_gcp=True))
    return [submit_j, copy_j]


def clean_ped_family_id(family_id: str) -> str:
    """
    Takes a family ID from the pedigree and cleans it up
    If the family ID already conforms to expectations, no action
    If the family ID fails, replace all non-alphanumeric/non-underscore characters with underscores

    >>> clean_ped_family_id('family1')
    'family1'
    >>> clean_ped_family_id('family-1-dirty')
    'family_1_dirty'

    Args:
        family_id (str): line from the pedigree file, unsplit

    Returns:
        the same line with a transformed family id
    """

    # if the family id is not valid, replace failing characters with underscores
    return (
        family_id
        if re.match(
            PED_FAMILY_ID_REGEX,
            family_id,
        )
        else re.sub(
            r'[^A-Za-z0-9_]',
            '_',
            family_id,
        )
    )


def make_combined_ped(cohort: targets.Cohort | targets.MultiCohort, combined_ped_path: Path) -> None:
    """
    Create cohort + ref panel PED.
    Concatenating all samples across all datasets with ref panel

    See #578 - there are restrictions on valid characters in PED file

    Args:
        cohort ():
        combined_ped_path ():

    Returns:
        None
    """

    # first get standard pedigree
    ped_dicts = [sequencing_group.pedigree.get_ped_dict() for sequencing_group in cohort.get_sequencing_groups()]

    if not ped_dicts:
        raise ValueError(f'No pedigree data found for {cohort.id}')

    conf_ped_path = get_references(['ped_file'])['ped_file']

    with (
        combined_ped_path.open('w') as out,
        to_path(conf_ped_path).open() as ref_ped,
    ):
        # layer of family ID cleaning
        for ped_dict in ped_dicts:
            out.write(
                '\t'.join(
                    [
                        clean_ped_family_id(ped_dict['Family.ID']),
                        ped_dict['Individual.ID'],
                        ped_dict['Father.ID'],
                        ped_dict['Mother.ID'],
                        ped_dict['Sex'],
                        ped_dict['Phenotype'],
                    ]
                )
                + '\n'
            )

        # The ref panel PED doesn't have any header, so can safely concatenate:
        out.write(ref_ped.read())


def queue_annotate_strvctvre_job(
    input_vcf,
    output_path: str,
    job_attrs: dict,
    name: str = 'AnnotateVcfWithStrvctvre',
) -> 'BashJob':
    """

    Args:
        input_vcf (ResourceFile): part of a resource group with the corresponding index
        output_path ():
        job_attrs (dict): job attributes
        name (str): name of the job

    Returns:
        The Strvctvre job
    """

    job_attrs = job_attrs or {}
    strv_job = hail_batch.get_batch().new_bash_job('StrVCTVRE', job_attrs | {'tool': 'strvctvre'})

    strv_job.image(config.config_retrieve(['images', 'strvctvre']))
    strv_job.storage(config.config_retrieve(['resource_overrides', name, 'storage'], '10Gi'))
    strv_job.memory(config.config_retrieve(['resource_overrides', name, 'memory'], '16Gi'))

    phylop = hail_batch.get_batch().read_input(config.config_retrieve(['references', 'strvctvre_phylop']))

    strv_job.declare_resource_group(output={'vcf.gz': '{root}.vcf.gz', 'vcf.gz.tbi': '{root}.vcf.gz.tbi'})

    # run strvctvre
    strv_job.command(f'python StrVCTVRE.py -i {input_vcf} -o {strv_job.output["vcf.gz"]} -f vcf -p {phylop}')
    strv_job.command(f'tabix {strv_job.output["vcf.gz"]}')

    hail_batch.get_batch().write_output(strv_job.output, str(output_path).replace('.vcf.gz', ''))
    return strv_job


def check_for_cohort_overlaps(multicohort: targets.MultiCohort):
    """
    Check for overlapping cohorts in a MultiCohort.
    GATK-SV does not tolerate overlapping cohorts, so we check for this here.
    This is called once per MultiCohort, and raises an Exception if any overlaps are found

    Args:
        multicohort (MultiCohort): the MultiCohort to check
    """
    # placeholder for errors
    errors: list[str] = []
    sgs_per_cohort: dict[str, set[str]] = {}
    # grab all SG IDs per cohort
    for cohort in multicohort.get_cohorts():
        # shouldn't be possible, but guard against to be sure
        if cohort.id in sgs_per_cohort:
            raise ValueError(f'Cohort {cohort.id} already exists in {sgs_per_cohort}')

        # collect the SG IDs for this cohort
        sgs_per_cohort[cohort.id] = set(cohort.get_sequencing_group_ids())

    # pairwise iteration over cohort IDs
    for id1, id2 in itertools.combinations(sgs_per_cohort, 2):
        # if the IDs are the same, skip. Again, shouldn't be possible, but guard against to be sure
        if id1 == id2:
            continue
        # if there are overlapping SGs, raise an error
        if overlap := sgs_per_cohort[id1] & sgs_per_cohort[id2]:
            errors.append(f'Overlapping cohorts {id1} and {id2} have overlapping SGs: {overlap}')
    # upon findings any errors, raise an Exception and die
    if errors:
        raise ValueError('\n'.join(errors))


def check_paths_exist(input_dict: dict[str, Any]):
    """
    Check that all paths in the input_dict exist
    """
    invalid_paths = defaultdict(list)
    for k in ['counts', 'PE_files', 'SD_files', 'SR_files'] + [f'{caller}_vcfs' for caller in SV_CALLERS]:
        for str_path in input_dict[k]:
            path = to_path(str_path)
            if not path.exists():
                sg_id = path.name.split('/')[-1].split('.')[0]
                invalid_paths[sg_id].append(str_path)
    if invalid_paths:
        error_str = '\n'.join([f'{k}: {", ".join(v)}' for k, v in invalid_paths.items()])
        raise FileNotFoundError(f'The following paths do not exist:\n{error_str}')


@cache
def query_for_spicy_vcf(dataset: str) -> str | None:
    """
    query for the most recent previous SpiceUpSVIDs VCF
    the SpiceUpSVIDs Stage involves overwriting the generic sequential variant IDs
    with meaningful Identifiers, so we can track the same variant across different callsets

    Args:
        dataset (str): project to query for

    Returns:
        str, the path to the latest Spicy VCF
        or None, if there are no Spicy VCFs
    """

    # hot swapping to a string we can freely modify
    query_dataset = dataset

    if config.config_retrieve(['workflow', 'access_level']) == 'test' and 'test' not in query_dataset:
        query_dataset += '-test'

    result = query(VCF_QUERY, variables={'dataset': query_dataset})
    spice_by_date: dict[str, str] = {}
    for analysis in result['project']['analyses']:
        if analysis['output'] and analysis['output'].endswith('fresh_ids.vcf.bgz'):
            spice_by_date[analysis['timestampCompleted']] = analysis['output']

    if not spice_by_date:
        return None

    # return the latest, determined by a sort on timestamp
    # 2023-10-10... > 2023-10-09..., so sort on strings
    return spice_by_date[sorted(spice_by_date)[-1]]


@functools.cache
def write_dataset_sg_ids(dataset: targets.Dataset) -> Path:
    """
    For a given dataset, write all its SGs to a file.
    Make this path specific to the dataset and run, so we can use it in multiple jobs

    Args:
        dataset ():

    Returns:

    """
    sgids_list_path = dataset.tmp_prefix() / workflow.get_workflow().output_version / 'sv-sgid-list.txt'
    if config.config_retrieve(['workflow', 'dry_run'], False):
        return sgids_list_path

    with sgids_list_path.open('w') as f:
        for sgid in dataset.get_sequencing_group_ids():
            f.write(f'{sgid}\n')

    return sgids_list_path
