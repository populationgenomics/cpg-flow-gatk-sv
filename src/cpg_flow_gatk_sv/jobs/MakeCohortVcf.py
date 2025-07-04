from typing import TYPE_CHECKING

from cpg_flow import targets
from cpg_flow_gatk_sv import utils
from cpg_utils import Path, config, to_path

if TYPE_CHECKING:
    from hailtop.batch.job import BashJob


def create_makecohortvcf_jobs(
    multicohort: targets.MultiCohort,
    pedigree_input: str,
    gatherbatchevidence_outputs: dict[str, Path],
    genotypebatch_outputs: dict[str, Path],
    filterbatch_outputs: dict[str, Path],
    outputs: dict[str, Path],
) -> list['BashJob']:
    fasta_file = config.config_retrieve(['workflow', 'ref_fasta'])

    # get the names of all contained cohorts
    all_batch_names = [cohort.id for cohort in multicohort.get_cohorts()]

    pesr_vcfs = [genotypebatch_outputs[cohort]['genotyped_pesr_vcf'] for cohort in all_batch_names]
    depth_vcfs = [genotypebatch_outputs[cohort]['genotyped_depth_vcf'] for cohort in all_batch_names]
    sr_pass = [genotypebatch_outputs[cohort]['sr_bothside_pass'] for cohort in all_batch_names]
    sr_fail = [genotypebatch_outputs[cohort]['sr_background_fail'] for cohort in all_batch_names]
    depth_depth_cutoff = [
        genotypebatch_outputs[cohort]['trained_genotype_depth_depth_sepcutoff'] for cohort in all_batch_names
    ]
    filter_batch_cutoffs = [filterbatch_outputs[cohort]['cutoffs'] for cohort in all_batch_names]
    bincov_files = [gatherbatchevidence_outputs[cohort]['merged_bincov'] for cohort in all_batch_names]
    disc_files = [gatherbatchevidence_outputs[cohort]['merged_PE'] for cohort in all_batch_names]
    median_cov_files = [gatherbatchevidence_outputs[cohort]['median_cov'] for cohort in all_batch_names]

    track_names = config.config_retrieve(['references', 'clustering_track_names'])
    track_bed_files = [
        config.config_retrieve(['references', 'clustering_track_sr']),
        config.config_retrieve(['references', 'clustering_track_sd']),
        config.config_retrieve(['references', 'clustering_track_rm']),
    ]

    input_dict = {
        # not explicit, but these VCFs require indices
        'HERVK_reference': config.config_retrieve(['references', 'hervk_reference']),
        'LINE1_reference': config.config_retrieve(['references', 'line1_reference']),
        'allosome_fai': config.config_retrieve(['references', 'allosome_file']),
        'batches': all_batch_names,
        'bin_exclude': config.config_retrieve(['references', 'bin_exclude']),
        'bincov_files': bincov_files,
        'chr_x': 'chrX',
        'chr_y': 'chrY',
        'clean_vcf1b_records_per_shard': 10000,
        'clean_vcf5_records_per_shard': 5000,
        'clustering_config_part1': config.config_retrieve(['references', 'clustering_config_part1']),
        'clustering_config_part2': config.config_retrieve(['references', 'clustering_config_part2']),
        'cohort_name': multicohort.name,
        'contig_list': config.config_retrieve(['references', 'primary_contigs_list']),
        'cytobands': config.config_retrieve(['references', 'cytobands']),
        'depth_gt_rd_sep_files': depth_depth_cutoff,
        # not explicit, but these VCFs require indices
        'depth_vcfs': depth_vcfs,
        'disc_files': disc_files,
        'gatk_docker': config.config_retrieve(['images', 'gatk_docker']),
        'linux_docker': config.config_retrieve(['images', 'linux_docker']),
        'max_shard_size_resolve': 500,
        'max_shards_per_chrom_clean_vcf_step1': 200,
        'median_coverage_files': median_cov_files,
        'mei_bed': config.config_retrieve(['references', 'mei_bed']),
        'min_records_per_shard_clean_vcf_step1': 5000,
        'min_sr_background_fail_batches': 0.5,
        'pe_exclude_list': config.config_retrieve(['references', 'pesr_exclude_list']),
        'ped_file': str(pedigree_input),
        # not explicit, but these VCFs require indices
        'pesr_vcfs': pesr_vcfs,
        'primary_contigs_list': config.config_retrieve(['references', 'primary_contigs_list']),
        'random_seed': 0,
        'raw_sr_background_fail_files': sr_fail,
        'raw_sr_bothside_pass_files': sr_pass,
        'reference_dict': str(to_path(fasta_file).with_suffix('.dict')),
        'reference_fasta': fasta_file,
        'reference_fasta_fai': f'{fasta_file}.fai',
        'rf_cutoff_files': filter_batch_cutoffs,
        'samples_per_clean_vcf_step2_shard': 100,
        'stratification_config_part1': config.config_retrieve(['references', 'stratification_config_part1']),
        'stratification_config_part2': config.config_retrieve(['references', 'stratification_config_part2']),
        'sv_base_mini_docker': config.config_retrieve(['images', 'sv_base_mini_docker']),
        'sv_pipeline_docker': config.config_retrieve(['images', 'sv_pipeline_docker']),
        'sv_pipeline_qc_docker': config.config_retrieve(['images', 'sv_pipeline_qc_docker']),
        'track_bed_files': track_bed_files,
        'track_names': track_names,
    }

    return utils.add_gatk_sv_jobs(
        dataset=multicohort.analysis_dataset,
        wfl_name='MakeCohortVcf',
        input_dict=input_dict,
        expected_out_dict=outputs,
        labels={'stage': 'makecohortvcf', config.AR_GUID_NAME: config.try_get_ar_guid()},
        job_size=utils.CromwellJobSizes.MEDIUM,
    )
