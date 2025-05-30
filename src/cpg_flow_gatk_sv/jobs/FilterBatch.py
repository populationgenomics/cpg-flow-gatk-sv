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
    """

    Args:
        pedigree_input ():
        cohort ():
        batch_metrics_outputs ():
        cluster_batch_outputs ():
        outputs ():

    Returns:

    """
    input_dict = {
        'batch': cohort.id,
        'ped_file': str(pedigree_input),
        'evidence_metrics': batch_metrics_outputs['metrics'],
        'evidence_metrics_common': batch_metrics_outputs['metrics_common'],
        'outlier_cutoff_nIQR': '6',
    }

    for caller in [*utils.SV_CALLERS, 'depth']:
        input_dict[f'{caller}_vcf'] = cluster_batch_outputs[f'clustered_{caller}_vcf']

    # add the images required for this step
    input_dict |= {
        key: config.config_retrieve(['images', key])
        for key in [
            'sv_pipeline_docker',
            'sv_base_mini_docker',
            'linux_docker',
        ]
    }

    input_dict |= utils.get_references(['primary_contigs_list'])

    # runs for over an hour
    return utils.add_gatk_sv_jobs(
        dataset=cohort.dataset,
        wfl_name='FilterBatch',
        input_dict=input_dict,
        expected_out_dict=outputs,
        labels={'stage': 'filterbatch', config.AR_GUID_NAME: config.try_get_ar_guid()},
        job_size=utils.CromwellJobSizes.MEDIUM,
    )
