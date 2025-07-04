from typing import TYPE_CHECKING

from cpg_flow import targets
from cpg_flow_gatk_sv import utils
from cpg_utils import Path, config, to_path

if TYPE_CHECKING:
    from hailtop.batch.job import BashJob


def create_joinrawcalls_jobs(
    multicohort: targets.MultiCohort,
    pedigree: str,
    clusterbatch_outputs: dict[str, Path],
    outputs: dict[str, Path],
) -> list['BashJob']:
    """

    Args:
        multicohort ():
        pedigree ():
        clusterbatch_outputs ():
        outputs ():

    Returns:

    """

    fasta_file = utils.get_fasta_string()
    input_dict = {
        'FormatVcfForGatk.formatter_args': '--fix-end',
        'prefix': multicohort.name,
        'ped_file': pedigree,
        'reference_fasta': fasta_file,
        'reference_fasta_fai': f'{fasta_file}.fai',
        'reference_dict': str(to_path(fasta_file).with_suffix('.dict')),
        'contig_list': config.config_retrieve(['references', 'primary_contigs_list']),
        'gatk_docker': config.config_retrieve(['images', 'gatk_docker']),
        'sv_base_mini_docker': config.config_retrieve(['images', 'sv_base_mini_docker']),
        'sv_pipeline_docker': config.config_retrieve(['images', 'sv_pipeline_docker']),
    }

    # get the names of all contained cohorts
    all_batch_names = [cohort.id for cohort in multicohort.get_cohorts()]

    for caller in [*utils.SV_CALLERS, 'depth']:
        input_dict[f'clustered_{caller}_vcfs'] = [
            clusterbatch_outputs[cohort][f'clustered_{caller}_vcf'] for cohort in all_batch_names
        ]
        input_dict[f'clustered_{caller}_vcf_indexes'] = [
            clusterbatch_outputs[cohort][f'clustered_{caller}_vcf_index'] for cohort in all_batch_names
        ]

    return utils.add_gatk_sv_jobs(
        dataset=multicohort.analysis_dataset,
        wfl_name='JoinRawCalls',
        input_dict=input_dict,
        expected_out_dict=outputs,
        labels={'stage': 'joinrawcalls', config.AR_GUID_NAME: config.try_get_ar_guid()},
    )
