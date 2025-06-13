from typing import TYPE_CHECKING

from cpg_flow import targets
from cpg_flow_gatk_sv import utils
from cpg_utils import Path, config, to_path

if TYPE_CHECKING:
    from hailtop.batch.job import BashJob


def create_cluster_batch_jobs(
    pedigree_input: str,
    cohort: targets.Cohort,
    batch_evidence_outputs: dict[str, Path],
    outputs: dict[str, Path],
) -> list['BashJob']:
    """

    Args:
        pedigree_input ():
        cohort ():
        batch_evidence_outputs ():
        outputs ():

    Returns:

    """

    fasta_file = utils.get_fasta_string()

    input_dict = {
        'batch': cohort.id,
        'del_bed': str(batch_evidence_outputs['merged_dels']),
        'dup_bed': str(batch_evidence_outputs['merged_dups']),
        'ped_file': pedigree_input,
        'depth_exclude_overlap_fraction': 0.5,
        'depth_interval_overlap': 0.8,
        'depth_clustering_algorithm': 'SINGLE_LINKAGE',
        'pesr_interval_overlap': 0.1,
        'pesr_breakend_window': 300,
        'pesr_clustering_algorithm': 'SINGLE_LINKAGE',
        'reference_fasta': fasta_file,
        'reference_fasta_fai': f'{fasta_file}.fai',
        'reference_dict': str(to_path(fasta_file).with_suffix('.dict')),
    }

    for caller in utils.SV_CALLERS:
        input_dict[f'{caller}_vcf_tar'] = str(batch_evidence_outputs[f'std_{caller}_vcf_tar'])

    # add the images required for this step
    input_dict |= {
        key: config.config_retrieve(['images', key])
        for key in [
            'linux_docker',
            'gatk_docker',
            'sv_base_mini_docker',
            'sv_pipeline_docker',
        ]
    }

    input_dict |= utils.get_references(
        [
            {'contig_list': 'primary_contigs_list'},
            {'depth_exclude_intervals': 'depth_exclude_list'},
            {'pesr_exclude_intervals': 'pesr_exclude_list'},
        ],
    )

    billing_labels = {'stage': 'clusterbatch', config.AR_GUID_NAME: config.try_get_ar_guid()}

    # runs for approx 1 hour
    return utils.add_gatk_sv_jobs(
        dataset=cohort.dataset,
        wfl_name='ClusterBatch',
        input_dict=input_dict,
        expected_out_dict=outputs,
        labels=billing_labels,
        job_size=utils.CromwellJobSizes.MEDIUM,
    )
