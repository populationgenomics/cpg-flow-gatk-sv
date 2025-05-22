import json
from typing import TYPE_CHECKING

import loguru
from cpg_flow import workflow
from cpg_utils import Path, config, hail_batch, to_path

from cpg_flow_gatk_sv.scripts import sample_batching

if TYPE_CHECKING:
    from hailtop.batch.job import BashJob


def create_sample_batches(
    qc_tables: list[Path],
    tmp_prefix: Path,
    output_json: str,
) -> 'BashJob':
    """
    pass all the QC tables from this run into a python script which will integrate that
    information into a set of sample batches for MultiSample stages

    Args:
        qc_tables ():
        tmp_prefix ():
        output_json (str): destination

    Returns:
        The BashJob
    """
    sequencing_groups = {
        sequencing_group.id: sequencing_group.meta
        for sequencing_group in workflow.get_multicohort().get_sequencing_groups()
    }
    if len(sequencing_groups) < config.config_retrieve(['workflow', 'min_batch_size'], 100):
        loguru.logger.error('Too few sequencing groups to form batches')
        raise RuntimeError('too few samples to create batches')

    # write them to a json file in tmp
    sgs_json_path = to_path(tmp_prefix / 'sgs_meta.json')
    with sgs_json_path.open('w') as f:
        json.dump(sequencing_groups, f)

    job = hail_batch.get_batch().new_bash_job('Create Sample Batches')
    job.image(config.config_retrieve(['workflow', 'driver_image']))
    job.command(f"""
        python3 {sample_batching.__file__} \
            --qc_tables {' '.join(qc_tables)} \
            --sgs_json {sgs_json_path} \
            --output_json {output_json}
    """)
    return job
