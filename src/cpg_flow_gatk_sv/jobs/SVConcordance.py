from typing import TYPE_CHECKING

from cpg_flow import targets
from cpg_flow_gatk_sv import utils
from cpg_utils import Path, config, to_path

if TYPE_CHECKING:
    from hailtop.batch.job import BashJob


def create_sv_concordance_jobs(
    multicohort: targets.MultiCohort,
    formatvcf_output: str,
    joinrawcalls_output: str,
    outputs: dict[str, Path],
) -> list['BashJob']:
    input_dict = {
        'output_prefix': multicohort.name,
        'reference_dict': str(to_path(config.config_retrieve(['workflow', 'ref_fasta'])).with_suffix('.dict')),
        'eval_vcf': formatvcf_output,
        'truth_vcf': joinrawcalls_output,
        'contig_list': config.config_retrieve(['references', 'primary_contigs_list']),
        'gatk_docker': config.config_retrieve(['images', 'gatk_docker']),
        'sv_base_mini_docker': config.config_retrieve(['images', 'sv_base_mini_docker']),
    }

    return utils.add_gatk_sv_jobs(
        dataset=multicohort.analysis_dataset,
        wfl_name='SVConcordance',
        input_dict=input_dict,
        expected_out_dict=outputs,
        labels={'stage': 'svconcordance', config.AR_GUID_NAME: config.try_get_ar_guid()},
    )
