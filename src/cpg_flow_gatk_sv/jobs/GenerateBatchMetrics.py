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
    """

    Args:
        pedigree_input ():
        cohort ():
        gather_batch_evidence_outputs ():
        cluster_batch_outputs ():
        outputs ():

    Returns:

    """

    fasta_file = utils.get_fasta_string()

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
        'ref_dict': str(to_path(fasta_file).with_suffix('.dict')),
    }

    for caller in [*utils.SV_CALLERS, 'depth']:
        input_dict[f'{caller}_vcf'] = cluster_batch_outputs[f'clustered_{caller}_vcf']

    # add the images required for this step
    input_dict |= {
        key: config.config_retrieve(['images', key])
        for key in [
            'sv_pipeline_docker',
            'sv_base_mini_docker',
            'sv_base_docker',
            'linux_docker',
        ]
    }

    input_dict |= utils.get_references(
        [
            'primary_contigs_list',
            'rmsk',
            'segdups',
            {'autosome_contigs': 'autosome_file'},
            {'allosome_contigs': 'allosome_file'},
        ],
    )

    # runs for approx 4-5 hours
    return utils.add_gatk_sv_jobs(
        dataset=cohort.dataset,
        wfl_name='GenerateBatchMetrics',
        input_dict=input_dict,
        expected_out_dict=outputs,
        labels={'stage': 'generatebatchmetrics', config.AR_GUID_NAME: config.try_get_ar_guid()},
        job_size=utils.CromwellJobSizes.MEDIUM,
    )
