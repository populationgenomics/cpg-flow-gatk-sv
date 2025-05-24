"""
Common methods for all GATK-SV workflows
"""

import re
from enum import Enum
from functools import cache
from os.path import join
from random import randint
from typing import TYPE_CHECKING, Any
import itertools

import loguru

from cpg_flow import targets, utils
from cpg_utils import Path, config, cromwell, hail_batch, to_path

if TYPE_CHECKING:
    from hailtop.batch.job import BashJob

GATK_SV_COMMIT = 'dc145a52f76a6f425ac3f481171040e78c0cfeea'
SV_CALLERS = ['manta', 'wham', 'scramble']

_FASTA_STRING = None
PED_FAMILY_ID_REGEX = re.compile(r'(^[A-Za-z0-9_]+$)')


class CromwellJobSizes(Enum):
    """
    Enum for polling intervals
    """

    SMALL = 'small'
    MEDIUM = 'medium'
    LARGE = 'large'


@cache
def get_sv_callers():
    if only_jobs := config.config_retrieve(['workflow', 'GatherSampleEvidence', 'only_jobs'], None):
        callers = [caller for caller in SV_CALLERS if caller in only_jobs]
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


def get_images(keys: list[str], allow_missing=False) -> dict[str, str]:
    """
    Dict of WDL inputs with docker image paths.

    Args:
        keys (list): all the images to get
        allow_missing (bool): if False, require all query keys to be found

    Returns:
        dict of image keys to image paths
        or AssertionError
    """
    image_keys = config.config_retrieve(['images']).keys()

    if not allow_missing:
        query_keys = set(keys)
        if not query_keys.issubset(image_keys):
            raise KeyError(f'Unknown image keys: {query_keys - image_keys}')

    return {k: config.image_path(k) for k in image_keys if k in keys}


def get_references(keys: list[str | dict[str, str]]) -> dict[str, str | list[str]]:
    """
    Dict of WDL inputs with reference file paths.
    """
    res: dict[str, str | list[str]] = {}
    for key in keys:
        # Keys can be maps (e.g. {'MakeCohortVcf.cytobands': 'cytoband'})
        ref_d_key = next(iter(key.values())) if isinstance(key, dict) else key

        # e.g. GATKSVPipelineBatch.rmsk -> rmsk
        ref_d_key = ref_d_key.split('.')[-1]
        try:
            res[key] = config.reference_path(f'gatk_sv/{ref_d_key}')  # type: ignore[index]
        except (KeyError, config.ConfigError):
            res[key] = config.reference_path(f'broad/{ref_d_key}')  # type: ignore[index]

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


def get_ref_panel(keys: list[str] | None = None) -> dict:
    # mandatory config entry
    ref_panel_samples = config.config_retrieve(['sv_ref_panel', 'ref_panel_samples'])
    return {
        k: v
        for k, v in {
            'ref_panel_samples': ref_panel_samples,
            'ref_panel_bincov_matrix': config.reference_path('broad/ref_panel_bincov_matrix'),
            'contig_ploidy_model_tar': config.reference_path('gatk_sv/contig_ploidy_model_tar'),
            'gcnv_model_tars': [
                str(config.reference_path('gatk_sv/model_tar_tmpl')).format(shard=i)
                for i in range(config.config_retrieve(['sv_ref_panel', 'model_tar_cnt']))
            ],
            'ref_panel_PE_files': [
                config.reference_path('gatk_sv/ref_panel_PE_file_tmpl').format(sample=s) for s in ref_panel_samples
            ],
            'ref_panel_SR_files': [
                config.reference_path('gatk_sv/ref_panel_SR_file_tmpl').format(sample=s) for s in ref_panel_samples
            ],
            'ref_panel_SD_files': [
                config.reference_path('gatk_sv/ref_panel_SD_file_tmpl').format(sample=s) for s in ref_panel_samples
            ],
        }.items()
        if not keys or k in keys
    }


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
            )

        # The ref panel PED doesn't have any header, so can safely concatenate:
        out.write(ref_ped.read())


def queue_annotate_sv_jobs(
    multicohort: targets.MultiCohort,
    prefix: Path,
    input_vcf: Path,
    outputs: dict,
    labels: dict[str, str] | None = None,
) -> list['BashJob']:
    """
    Helper function to queue jobs for SV annotation
    Enables common access to the same Annotation WDL for CNV & SV
    """
    input_dict: dict[str, Any] = {
        'vcf': input_vcf,
        'prefix': multicohort.name,
        'ped_file': make_combined_ped(multicohort, prefix),
        'sv_per_shard': 5000,
        'external_af_population': config.config_retrieve(['references', 'gatk_sv', 'external_af_population']),
        'external_af_ref_prefix': config.config_retrieve(['references', 'gatk_sv', 'external_af_ref_bed_prefix']),
        'external_af_ref_bed': config.config_retrieve(['references', 'gnomad_sv']),
        'use_hail': False,
    }

    input_dict |= get_references(
        [
            'noncoding_bed',
            'protein_coding_gtf',
            {'contig_list': 'primary_contigs_list'},
        ],
    )

    # images!
    input_dict |= get_images(['sv_pipeline_docker', 'sv_base_mini_docker', 'gatk_docker'])
    return add_gatk_sv_jobs(
        dataset=multicohort.analysis_dataset,
        wfl_name='AnnotateVcf',
        input_dict=input_dict,
        expected_out_dict=outputs,
        labels=labels,
    )


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

    strv_job.image(config.image_path('strvctvre'))
    strv_job.storage(config.config_retrieve(['resource_overrides', name, 'storage'], '10Gi'))
    strv_job.memory(config.config_retrieve(['resource_overrides', name, 'memory'], '16Gi'))

    strvctvre_phylop = get_references(['strvctvre_phylop'])['strvctvre_phylop']

    local_phylop = hail_batch.get_batch().read_input(strvctvre_phylop)

    strv_job.declare_resource_group(output={'vcf.gz': '{root}.vcf.gz', 'vcf.gz.tbi': '{root}.vcf.gz.tbi'})

    # run strvctvre
    strv_job.command(
        f'python StrVCTVRE.py -i {input_vcf} -o {strv_job.output["vcf.gz"]} -f vcf -p {local_phylop}',
    )
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
