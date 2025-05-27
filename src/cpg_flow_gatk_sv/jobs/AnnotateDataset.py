from typing import TYPE_CHECKING

from cpg_utils import config, Path, hail_batch

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
    """
    parser = argparse.ArgumentParser(description='Annotate SV dataset MatrixTable')
    parser.add_argument('--mt', required=True, help='Path to the input MatrixTable')
    parser.add_argument('--dataset_mt', required=True, help='Path to write the single-dataset MT')
    parser.add_argument('--out_mt', required=True, help='Path to write the final MatrixTable')
    parser.add_argument(
        '--sample_id_file',
        required=True,
        help='Path to a file containing sample IDs to include in the subset, one per line',
    )
    parser.add_argument(
        '--exclusion_file',
        required=False,
        help='Path to a file containing sample IDs to exclude from the subset, one per line',
    )
    args = parser.parse_args()
    """
    job = hail_batch.get_batch().new_bash_job('Annotate dataset', job_attrs)
    job.image(config.config_retrieve(['workflow', 'driver_image']))
    job.command(f"""
    python3 -m cpg_flow_gatk_sv.scripts.annotate_dataset \
            --mt {mt} \
            --dataset_mt {dataset_mt} \
            --out_mt {mt_out} \
            --sample_id_file {sgid_file} \
            --exclusion_file {exclusion_file}
    """)
    return job
