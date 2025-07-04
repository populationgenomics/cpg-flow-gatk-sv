from typing import TYPE_CHECKING

from cpg_flow import targets
from cpg_flow_gatk_sv import utils
from cpg_utils import Path, config, to_path

if TYPE_CHECKING:
    from hailtop.batch.job import BashJob


def create_generate_batch_metrics_jobs(
    pedigree_input: str,
    cohort: targets.Cohort,
    gather_batch_evidence_outputs: dict[str, Path],
    cluster_batch_outputs: dict[str, Path],
    outputs: dict[str, Path],
) -> list['BashJob']:
    input_dict = {
        'batch': cohort.id,
        'baf_metrics': gather_batch_evidence_outputs['merged_BAF'],
        'discfile': gather_batch_evidence_outputs['merged_PE'],
        'coveragefile': gather_batch_evidence_outputs['merged_bincov'],
        'splitfile': gather_batch_evidence_outputs['merged_SR'],
        'medianfile': gather_batch_evidence_outputs['median_cov'],
        'BAF_split_size': 10000,
        'RD_split_size': 10000,
        'PE_split_size': 10000,
        'SR_split_size': 1000,
        'common_cnv_size_cutoff': 5000,
        'ped_file': pedigree_input,
        'ref_dict': str(to_path(config.config_retrieve(['workflow', 'ref_fasta'])).with_suffix('.dict')),
        'primary_contigs_list': config.config_retrieve(['references', 'primary_contigs_list']),
        'rmsk': config.config_retrieve(['references', 'rmsk']),
        'segdups': config.config_retrieve(['references', 'segdups']),
        'autosome_contigs': config.config_retrieve(['references', 'autosome_file']),
        'allosome_contigs': config.config_retrieve(['references', 'allosome_file']),
        'sv_pipeline_docker': config.config_retrieve(['images', 'sv_pipeline_docker']),
        'sv_base_mini_docker': config.config_retrieve(['images', 'sv_base_mini_docker']),
        'sv_base_docker': config.config_retrieve(['images', 'sv_base_docker']),
        'linux_docker': config.config_retrieve(['images', 'linux_docker']),
    }

    for caller in [*utils.SV_CALLERS, 'depth']:
        input_dict[f'{caller}_vcf'] = cluster_batch_outputs[f'clustered_{caller}_vcf']

    # runs for approx 4-5 hours
    return utils.add_gatk_sv_jobs(
        dataset=cohort.dataset,
        wfl_name='GenerateBatchMetrics',
        input_dict=input_dict,
        expected_out_dict=outputs,
        labels={'stage': 'generatebatchmetrics', config.AR_GUID_NAME: config.try_get_ar_guid()},
        job_size=utils.CromwellJobSizes.MEDIUM,
    )
