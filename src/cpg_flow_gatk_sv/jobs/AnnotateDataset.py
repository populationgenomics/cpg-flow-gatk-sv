from typing import TYPE_CHECKING

from cpg_utils import Path, config, hail_batch

if TYPE_CHECKING:
    from hailtop.batch.job import BashJob


def create_annotate_dataset_jobs(
    mt: str,
    sgid_file: Path,
    mt_out: str,
    dataset_mt: Path,
    exclusion_file: Path,
    job_attrs: dict,
) -> 'BashJob':
    job = hail_batch.get_batch().new_bash_job('Annotate dataset', job_attrs)
    job.image(config.config_retrieve(['workflow', 'driver_image']))
    job.command(f"""
    python3 -m cpg_flow_gatk_sv.scripts.annotate_dataset \
            --mt {mt} \
            --dataset_mt {dataset_mt} \
            --output {mt_out} \
            --sample_id_file {sgid_file} \
            --exclusion_file {exclusion_file}
    """)
    return job
