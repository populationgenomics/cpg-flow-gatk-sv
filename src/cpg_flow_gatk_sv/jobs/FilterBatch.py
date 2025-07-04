from typing import TYPE_CHECKING

from cpg_flow import targets
from cpg_flow_gatk_sv import utils
from cpg_utils import Path, config

if TYPE_CHECKING:
    from hailtop.batch.job import BashJob


def create_filterbatch_jobs(
    pedigree_input: str,
    cohort: targets.Cohort,
    batch_metrics_outputs: dict[str, Path],
    cluster_batch_outputs: dict[str, Path],
    outputs: dict[str, Path],
) -> list['BashJob']:
    input_dict = {
        'batch': cohort.id,
        'ped_file': str(pedigree_input),
        'evidence_metrics': batch_metrics_outputs['metrics'],
        'evidence_metrics_common': batch_metrics_outputs['metrics_common'],
        'outlier_cutoff_nIQR': '6',
        'linux_docker': config.config_retrieve(['images', 'linux_docker']),
        'sv_base_mini_docker': config.config_retrieve(['images', 'sv_base_mini_docker']),
        'sv_pipeline_docker': config.config_retrieve(['images', 'sv_pipeline_docker']),
        'primary_contigs_list': config.config_retrieve(['references', 'primary_contigs_list']),
    }

    for caller in [*utils.SV_CALLERS, 'depth']:
        input_dict[f'{caller}_vcf'] = cluster_batch_outputs[f'clustered_{caller}_vcf']

    # runs for over an hour
    return utils.add_gatk_sv_jobs(
        dataset=cohort.dataset,
        wfl_name='FilterBatch',
        input_dict=input_dict,
        expected_out_dict=outputs,
        labels={'stage': 'filterbatch', config.AR_GUID_NAME: config.try_get_ar_guid()},
        job_size=utils.CromwellJobSizes.MEDIUM,
    )
