from typing import TYPE_CHECKING

from cpg_utils import config, hail_batch

if TYPE_CHECKING:
    from hailtop.batch.job import BashJob


def create_combine_exclusion_lists_job(
    file_list: list[str],
    output: str,
) -> 'BashJob':
    job = hail_batch.get_batch().new_bash_job('Concatenate all sample exclusion files')
    job.image(config.config_retrieve(['workflow', 'driver_image']))
    hail_batch.authenticate_cloud_credentials_in_job(job)
    job.command(f'gcloud storage objects compose {" ".join(file_list)} {output}')

    return job
