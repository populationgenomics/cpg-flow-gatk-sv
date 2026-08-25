from cpg_utils import Path, config, hail_batch


def cohort_to_vcf_job(
    dataset: str,
    input_mt: str,
    dataset_sgids: Path,
    exclusions: str,
    output_vcf: Path,
    job_attrs: dict,
):
    """Take the single-dataset MT, and write to a VCF."""

    job = hail_batch.get_batch().new_bash_job('VCF from dataset MT', job_attrs | {'tool': 'hail query'})
    job.image(config.config_retrieve(['workflow', 'driver_image']))
    job.command(f'python -m cpg_flow_gatk_sv.scripts.vcf_from_mt_subset --input {input_mt} --output {output_vcf!s}')

    job.command(f"""
    python3 -m cpg_flow_gatk_sv.scripts.register_with_exclusions \\
        --output {output_vcf!s} \\
        --dataset {dataset} \\
        --stage AnnotatedDatasetMtToSvVcf \\
        --atype custom \\
        --sgs {dataset_sgids!s} \\
        --exclusions {exclusions!s}
    """)
    return job
