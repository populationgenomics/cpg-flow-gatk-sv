from typing import TYPE_CHECKING

from cpg_flow import targets
from cpg_flow_gatk_sv import utils
from cpg_utils import Path, config

if TYPE_CHECKING:
    from hailtop.batch.job import BashJob


def create_filtergenotypes_jobs(
    multicohort: targets.MultiCohort,
    pedigree: str,
    ploidy_table: str,
    sv_conc_vcf: str,
    outputs: dict[str, Path],
) -> list['BashJob']:
    input_dict = {
        'output_prefix': multicohort.name,
        'vcf': sv_conc_vcf,
        'ploidy_table': ploidy_table,
        'ped_file': pedigree,
        'fmax_beta': config.config_retrieve(['references', 'fmax_beta'], 0.4),
        'recalibrate_gq_args': config.config_retrieve(['references', 'recalibrate_gq_args']),
        'sl_filter_args': config.config_retrieve(['references', 'sl_filter_args']),
        'primary_contigs_fai': config.config_retrieve(['references', 'primary_contigs_fai']),
        'gq_recalibrator_model_file': config.config_retrieve(['references', 'aou_filtering_model']),
        'gatk_docker': config.config_retrieve(['images', 'gq_recalibrator_docker']),
        'linux_docker': config.config_retrieve(['images', 'linux_docker']),
        'sv_base_mini_docker': config.config_retrieve(['images', 'sv_base_mini_docker']),
        'sv_pipeline_docker': config.config_retrieve(['images', 'sv_pipeline_docker']),
        'genome_tracks': config.config_retrieve(['references', 'genome_tracks']),
    }

    return utils.add_gatk_sv_jobs(
        dataset=multicohort.analysis_dataset,
        wfl_name='FilterGenotypes',
        input_dict=input_dict,
        expected_out_dict=outputs,
        labels={'stage': 'filtergenotypes', config.AR_GUID_NAME: config.try_get_ar_guid()},
    )
