from typing import TYPE_CHECKING

from cpg_utils import Path, config, hail_batch

if TYPE_CHECKING:
    from hailtop.batch.job import BashJob


def create_combine_exclusion_lists_job(
    file_list: list[str],
    output: Path,
) -> 'BashJob':
    job = hail_batch.get_batch().new_bash_job('Concatenate all sample exclusion files')
    job.image(config.config_retrieve(['workflow', 'driver_image']))
    hail_batch.authenticate_cloud_credentials_in_job(job)

    # Gcloud compose allows concatenation of up to 32 files, more than this needs chunking
    intermediate_files = []
    for i in range(0, len(file_list), 32):
        chunk_files = file_list[i : i + 32]
        intermediate_output = f'{output!s}_intermediate_{i // 32}'
        job.command(f'gcloud storage objects compose {" ".join(chunk_files)} {intermediate_output}')
        intermediate_files.append(intermediate_output)

    job.command(f'gcloud storage objects compose {" ".join(intermediate_files)} {output!s}')

    return job
