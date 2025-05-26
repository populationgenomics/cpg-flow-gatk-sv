from typing import TYPE_CHECKING

from cpg_flow import targets
from cpg_flow_gatk_sv import utils
from cpg_utils import Path, config, to_path

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
        'fmax_beta': config.config_retrieve(['references', 'gatk_sv', 'fmax_beta'], 0.4),
        'recalibrate_gq_args': config.config_retrieve(['references', 'gatk_sv', 'recalibrate_gq_args']),
        'sl_filter_args': config.config_retrieve(['references', 'gatk_sv', 'sl_filter_args']),
    }
    input_dict |= {
        key: config.config_retrieve(['images', key])
        for key in [
            'linux_docker',
            'sv_base_mini_docker',
            'sv_pipeline_docker',
        ]
    }

    # use a non-standard GATK image containing required filtering tool
    input_dict['gatk_docker'] = config.config_retrieve(['images', 'gq_recalibrator_docker'])
    input_dict |= utils.get_references(
        [
            {'gq_recalibrator_model_file': 'aou_filtering_model'},
            'primary_contigs_fai',
        ],
    )

    # something a little trickier - we need to get various genome tracks
    input_dict['genome_tracks'] = list(
        utils.get_references(config.config_retrieve(['references', 'gatk_sv', 'genome_tracks'], [])).values(),
    )

    return utils.add_gatk_sv_jobs(
        dataset=multicohort.analysis_dataset,
        wfl_name='FilterGenotypes',
        input_dict=input_dict,
        expected_out_dict=outputs,
        labels={'stage': 'filtergenotypes', config.AR_GUID_NAME: config.try_get_ar_guid()},
    )
