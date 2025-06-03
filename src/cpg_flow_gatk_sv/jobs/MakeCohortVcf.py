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
    fasta_file = utils.get_fasta_string()

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

    track_names = config.config_retrieve(['references', 'gatk_sv', 'clustering_track_names'])
    track_bed_files = [
        config.config_retrieve(['references', 'gatk_sv', 'clustering_track_sr']),
        config.config_retrieve(['references', 'gatk_sv', 'clustering_track_sd']),
        config.config_retrieve(['references', 'gatk_sv', 'clustering_track_rm']),
    ]

    input_dict = {
        'cohort_name': multicohort.name,
        'batches': all_batch_names,
        'ped_file': str(pedigree_input),
        'reference_fasta': fasta_file,
        'reference_fasta_fai': f'{fasta_file}.fai',
        'reference_dict': str(to_path(fasta_file).with_suffix('.dict')),
        'chr_x': 'chrX',
        'chr_y': 'chrY',
        'min_sr_background_fail_batches': 0.5,
        'max_shard_size_resolve': 500,
        'max_shards_per_chrom_clean_vcf_step1': 200,
        'min_records_per_shard_clean_vcf_step1': 5000,
        'clean_vcf1b_records_per_shard': 10000,
        'samples_per_clean_vcf_step2_shard': 100,
        'clean_vcf5_records_per_shard': 5000,
        'random_seed': 0,
        # not explicit, but these VCFs require indices
        'pesr_vcfs': pesr_vcfs,
        'depth_vcfs': depth_vcfs,
        'disc_files': disc_files,
        'bincov_files': bincov_files,
        'raw_sr_bothside_pass_files': sr_pass,
        'raw_sr_background_fail_files': sr_fail,
        'depth_gt_rd_sep_files': depth_depth_cutoff,
        'median_coverage_files': median_cov_files,
        'rf_cutoff_files': filter_batch_cutoffs,
        'track_names': track_names,
        'track_bed_files': track_bed_files,
    }

    input_dict |= utils.get_references(
        [
            'bin_exclude',
            'clustering_config_part1',
            'clustering_config_part2',
            'mei_bed',
            'stratification_config_part1',
            'stratification_config_part2',
            # same attr, two names
            'primary_contigs_list',
            {'allosome_fai': 'allosome_file'},
            {'contig_list': 'primary_contigs_list'},
            {'cytobands': 'cytoband'},
            {'HERVK_reference': 'hervk_reference'},
            {'LINE1_reference': 'line1_reference'},
            {'pe_exclude_list': 'pesr_exclude_list'},
        ],
    )

    # add the images required for this step
    input_dict |= {
        key: config.config_retrieve(['images', key])
        for key in [
            'gatk_docker',
            'sv_pipeline_docker',
            'sv_pipeline_qc_docker',
            'sv_base_mini_docker',
            'linux_docker',
        ]
    }

    return utils.add_gatk_sv_jobs(
        dataset=multicohort.analysis_dataset,
        wfl_name='MakeCohortVcf',
        input_dict=input_dict,
        expected_out_dict=outputs,
        labels={'stage': 'makecohortvcf', config.AR_GUID_NAME: config.try_get_ar_guid()},
        job_size=utils.CromwellJobSizes.MEDIUM,
    )
