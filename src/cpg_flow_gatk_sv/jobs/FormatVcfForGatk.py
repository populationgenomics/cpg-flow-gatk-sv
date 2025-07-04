from typing import TYPE_CHECKING

from cpg_flow import targets
from cpg_flow_gatk_sv import utils
from cpg_utils import Path, config

if TYPE_CHECKING:
    from hailtop.batch.job import BashJob


def create_formatvcf_jobs(
    multicohort: targets.MultiCohort,
    pedigree_input: str,
    vcf_input: str,
    outputs: dict[str, Path],
) -> list['BashJob']:
    input_dict = {
        'prefix': multicohort.name,
        'vcf': vcf_input,
        'ped_file': pedigree_input,
        'contig_list': config.config_retrieve(['references', 'primary_contigs_list']),
        'sv_pipeline_docker': config.config_retrieve(['images', 'sv_pipeline_docker']),
        'sv_base_mini_docker': config.config_retrieve(['images', 'sv_base_mini_docker']),
    }

    return utils.add_gatk_sv_jobs(
        dataset=multicohort.analysis_dataset,
        wfl_name='FormatVcfForGatk',
        input_dict=input_dict,
        expected_out_dict=outputs,
        labels={'stage': 'formatvcfforgatk', config.AR_GUID_NAME: config.try_get_ar_guid()},
        job_size=utils.CromwellJobSizes.SMALL,
    )
