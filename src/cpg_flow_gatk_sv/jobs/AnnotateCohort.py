from typing import TYPE_CHECKING

from cpg_utils import Path, config, hail_batch

if TYPE_CHECKING:
    from hailtop.batch.job import BashJob


def create_annotate_cohort_job(
    vcf: Path,
    checkpoint: Path,
    out_mt: Path,
    job_attrs: dict,
) -> 'BashJob':
    job = hail_batch.get_batch().new_bash_job('Annotate cohort', job_attrs)
    job.image(config.config_retrieve(['workflow', 'driver_image']))
    gencode_gtf_local = hail_batch.get_batch().read_input(config.config_retrieve(['workflow', 'gencode_gtf_file']))
    job.storage('10Gi')
    job.command(f"""
        python -m cpg_flow_gatk_sv.scripts.annotate_cohort \
        --vcf {vcf!s}
        --output {out_mt!s}
        --gencode_gtf {gencode_gtf_local}
        --checkpoint {checkpoint!s}
    """)
    return job
