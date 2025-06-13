import functools
from typing import TYPE_CHECKING

import loguru
from google.api_core import exceptions

from cpg_flow import targets
from cpg_utils import Path, cloud, config, hail_batch

if TYPE_CHECKING:
    from hailtop.batch.job import BashJob


@functools.cache
def es_password() -> str:
    """
    Get Elasticsearch password. Moved into a separate method to simplify
    mocking in tests.
    """
    return cloud.read_secret(
        project_id=config.config_retrieve(['elasticsearch', 'password_project_id']),
        secret_name=config.config_retrieve(['elasticsearch', 'password_secret_id']),
        fail_gracefully=False,
    )


def create_mt_to_es_job(
    dataset: targets.Dataset,
    mt_path: str,
    outputs: dict[str, Path],
    job_attrs: dict,
) -> 'BashJob | None':
    """
    Create a job to export a MatrixTable to Elasticsearch.
    """

    # try to generate a password here - we'll find out inside the script anyway, but
    # by that point we'd already have localised the MT, wasting time and money
    try:
        _es_password_string = es_password()
    except exceptions.PermissionDenied:
        loguru.logger.warning(f'No permission to access ES password, skipping for {dataset.name}')
        return None
    except KeyError:
        loguru.logger.warning(f'ES section not in config, skipping for {dataset.name}')
        return None

    job = hail_batch.get_batch().new_bash_job(f'Export MT to ES: {dataset.name}', job_attrs)
    job.image(config.config_retrieve(['workflow', 'driver_image'])).cpu(4).memory('lowmem').storage('10Gi')

    # and just the name, used after localisation
    mt_name = mt_path.split('/')[-1]

    # localise the MT
    job.command(f'gcloud --no-user-output-enabled storage cp -r {mt_path} $BATCH_TMPDIR')

    job.command(f"""
        python3 -m cpg_flow_gatk_sv.scripts.mt_to_es_export \
            --mt "${{BATCH_TMPDIR}}/{mt_name}" \
            --es_index {outputs['index_name']!s} \
            --flag {outputs['flag']!s} \
    """)

    return job
