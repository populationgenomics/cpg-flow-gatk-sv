from typing import TYPE_CHECKING

from cpg_flow import targets
from cpg_flow_gatk_sv import utils
from cpg_utils import Path, config, to_path

if TYPE_CHECKING:
    from hailtop.batch.job import BashJob


def submit_gatherbatchevidence_jobs(
    pedigree_input: Path,
    cohort: targets.Cohort,
    train_gcnv: dict[str, Path | list[Path]],
    outputs: dict[str, Path],
) -> list['BashJob']:
    """

    Args:
        pedigree_input ():
        cohort ():
        train_gcnv (dict[str, Path | list[Path]]): output from TrainGCNV stage, containing the model tar files
        outputs (dict[str, Path]): expected outputs of this stage

    Returns:

    """

    sequencing_groups = cohort.get_sequencing_groups(only_active=True)

    fasta_file = utils.get_fasta_string()

    input_dict = {
        'batch': cohort.id,
        'samples': [sg.id for sg in sequencing_groups],
        'ped_file': str(pedigree_input),
        'counts': [
            str(sequencing_group.make_sv_evidence_path / f'{sequencing_group.id}.coverage_counts.tsv.gz')
            for sequencing_group in sequencing_groups
        ],
        'SR_files': [
            str(sequencing_group.make_sv_evidence_path / f'{sequencing_group.id}.sr.txt.gz')
            for sequencing_group in sequencing_groups
        ],
        'PE_files': [
            str(sequencing_group.make_sv_evidence_path / f'{sequencing_group.id}.pe.txt.gz')
            for sequencing_group in sequencing_groups
        ],
        'SD_files': [
            str(sequencing_group.make_sv_evidence_path / f'{sequencing_group.id}.sd.txt.gz')
            for sequencing_group in sequencing_groups
        ],
        'ref_copy_number_autosomal_contigs': 2,
        'allosomal_contigs': ['chrX', 'chrY'],
        'gcnv_qs_cutoff': 30,
        'min_svsize': 50,
        'run_matrix_qc': True,
        'matrix_qc_distance': 1000000,
        'ref_dict': str(to_path(fasta_file).with_suffix('.dict')),
        'genome_file': config.config_retrieve(['references', 'genome_file']),
        'primary_contigs_fai': config.config_retrieve(['references', 'primary_contigs_fai']),
        'sd_locs_vcf': config.config_retrieve(['references', 'dbsnp_vcf']),
        'cnmops_chrom_file': config.config_retrieve(['references', 'autosome_file']),
        'cnmops_exclude_list': config.config_retrieve(['references', 'cnmops_exclude_list']),
        'cnmops_allo_file': config.config_retrieve(['references', 'allosome_file']),
        'cytoband': config.config_retrieve(['references', 'cytobands']),
        'mei_bed': config.config_retrieve(['references', 'mei_bed']),
        'sv_base_docker': config.config_retrieve(['images', 'sv_base_docker']),
        'sv_base_mini_docker': config.config_retrieve(['images', 'sv_base_mini_docker']),
        'sv_pipeline_docker': config.config_retrieve(['images', 'sv_pipeline_docker']),
        'sv_pipeline_qc_docker': config.config_retrieve(['images', 'sv_pipeline_qc_docker']),
        'linux_docker': config.config_retrieve(['images', 'linux_docker']),
        'condense_counts_docker': config.config_retrieve(['images', 'condense_counts_docker']),
        'gatk_docker': config.config_retrieve(['images', 'gatk_docker']),
        'cnmops_docker': config.config_retrieve(['images', 'cnmops_docker']),
    }

    for caller in utils.SV_CALLERS:
        input_dict[f'{caller}_vcfs'] = [
            str(sequencing_group.make_sv_evidence_path / f'{sequencing_group.id}.{caller}.vcf.gz')
            for sequencing_group in sequencing_groups
        ]

    # The inputs are from a different workflow, so we need to check if they exist manually
    # do this before adding the TrainGCNV outputs, which are generated in this workflow, so may not exist yet
    utils.check_paths_exist(input_dict)

    # reference panel gCNV models - a mix of canonical and cohort-specific input data
    ref_panel_samples = config.config_retrieve(['sv_ref_panel', 'ref_panel_samples'])
    input_dict |= {
        'ref_panel_samples': ref_panel_samples,
        'ref_panel_bincov_matrix': config.reference_path('broad/ref_panel_bincov_matrix'),
        'contig_ploidy_model_tar': str(train_gcnv['cohort_contig_ploidy_model_tar']),
        'gcnv_model_tars': [str(x) for x in train_gcnv['cohort_gcnv_model_tars']],
        'ref_panel_PE_files': [
            config.config_retrieve(['references', 'ref_panel_PE_file_tmpl']).format(sample=s) for s in ref_panel_samples
        ],
        'ref_panel_SR_files': [
            config.config_retrieve(['references', 'ref_panel_SR_file_tmpl']).format(sample=s) for s in ref_panel_samples
        ],
        'ref_panel_SD_files': [
            config.config_retrieve(['references', 'ref_panel_SD_file_tmpl']).format(sample=s) for s in ref_panel_samples
        ],
    }

    # this step runs for approximately 15 hours
    return utils.add_gatk_sv_jobs(
        dataset=cohort.dataset,
        wfl_name='GatherBatchEvidence',
        input_dict=input_dict,
        expected_out_dict=outputs,
        labels={'stage': 'gatherbatchevidence', config.AR_GUID_NAME: config.try_get_ar_guid()},
        job_size=utils.CromwellJobSizes.LARGE,
    )
