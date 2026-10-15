from typing import TYPE_CHECKING

from cpg_flow import targets
from cpg_utils import Path, config, hail_batch

if TYPE_CHECKING:
    from hailtop.batch.job import BashJob


def write_sample_map(dataset: targets.Dataset, sg_to_participant: dict[str, str], sample_map_file: Path) -> None:
    """Write a TSV of SG ID, participant external ID and family ID, for every SG to export."""
    if config.config_retrieve(['workflow', 'dry_run'], False):
        return

    family_by_sg = {sg.id: sg.pedigree.fam_id for sg in dataset.get_sequencing_groups()}
    with sample_map_file.open('w') as handle:
        for sg_id, participant_id in sorted(sg_to_participant.items()):
            handle.write(f'{sg_id}\t{participant_id}\t{family_by_sg[sg_id]}\n')


def create_seqr_sv_vcf_jobs(
    dataset: targets.Dataset,
    input_vcf: str,
    sg_to_participant: dict[str, str],
    sample_map_file: Path,
    exclusions: str,
    output_vcf: Path,
    job_attrs: dict[str, str],
) -> list['BashJob']:
    """
    Subset the multicohort SV VCF to this dataset, rewrite it into a format seqr can load, then bgzip and index

    Args:
        dataset (Dataset): the dataset being exported
        input_vcf (str): the multicohort SpiceUpSvIds VCF
        sg_to_participant (dict): {SG ID: participant external ID} for every SG to export
        sample_map_file (Path): where to write the SG ID -> participant ID -> family ID TSV
        exclusions (str): the combined list of SG IDs removed by QC
        output_vcf (Path): where to write the seqr-ready VCF
        job_attrs (dict): attributes for the jobs
    """
    write_sample_map(dataset, sg_to_participant, sample_map_file)

    batch = hail_batch.get_batch()
    local_vcf = batch.read_input(input_vcf)
    local_sample_map = batch.read_input(sample_map_file)
    local_exclusions = batch.read_input(exclusions)

    # drop other datasets' samples with bcftools first, as the multicohort VCF is very large uncompressed
    # --force-samples allows for excluded SGs, which aren't in the VCF. The rewrite fails on any other missing SG
    # --no-update keeps the cohort-level AC/AN, consistent with AF, N_HET and N_HOMALT
    subset_job = batch.new_bash_job(f'SeqrSvVcf subset: {dataset.name}', job_attrs | {'tool': 'bcftools'})
    subset_job.image(config.config_retrieve(['images', 'bcftools']))
    subset_job.cpu(1).memory('highmem').storage('10Gi')
    subset_job.command(
        f"""
        cut -f1 {local_sample_map} > sg_ids.txt
        bcftools view \\
            --force-samples \\
            --no-update \\
            -S sg_ids.txt \\
            -Oz \\
            -o {subset_job.output} \\
            {local_vcf}
        """
    )

    gnomad_field = f'{config.config_retrieve(["references", "external_af_ref_bed_prefix"])}_SVID'
    rewrite_job = batch.new_bash_job(f'SeqrSvVcf rewrite: {dataset.name}', job_attrs)
    rewrite_job.image(config.config_retrieve(['workflow', 'driver_image']))
    rewrite_job.storage('10Gi')
    rewrite_job.command(
        f"""
        python -m cpg_flow_gatk_sv.scripts.seqr_sv_vcf \\
            --input {subset_job.output} \\
            --sample_map {local_sample_map} \\
            --exclusions {local_exclusions} \\
            --gnomad_field {gnomad_field} \\
            --output {rewrite_job.output}
        """
    )

    compress_job = batch.new_bash_job(f'SeqrSvVcf bgzip and tabix: {dataset.name}', job_attrs | {'tool': 'bcftools'})
    compress_job.image(config.config_retrieve(['images', 'bcftools']))
    compress_job.storage('10Gi')
    compress_job.declare_resource_group(
        output={
            'vcf.bgz': '{root}.vcf.bgz',
            'vcf.bgz.tbi': '{root}.vcf.bgz.tbi',
        }
    )
    compress_job.command(f'bcftools view {rewrite_job.output} -Oz -W=tbi -o {compress_job.output["vcf.bgz"]}')

    batch.write_output(compress_job.output, str(output_vcf).removesuffix('.vcf.bgz'))

    return [subset_job, rewrite_job, compress_job]
