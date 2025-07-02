from typing import TYPE_CHECKING

from cpg_flow_gatk_sv import utils
from cpg_utils import config, hail_batch

if TYPE_CHECKING:
    from hailtop.batch.job import BashJob


def create_generate_ploidy_jobs(
    pedigree: str,
    output: str,
) -> 'BashJob':
    """

    Args:
        pedigree ():
        output ():

    Returns:

    """

    job = hail_batch.get_batch().new_bash_job('Run create_ploidy_table script')
    job.image(config.config_retrieve(['workflow', 'driver_image']))

    contig_path = utils.get_references(['primary_contigs_list'])['primary_contigs_list']

    job.command(f"""
    python3 -m cpg_flow_gatk_sv.scripts.ploidy_table_from_ped \
            --ped {pedigree} \
            --output {job.output} \
            --contigs {contig_path}
    """)

    hail_batch.get_batch().write_output(job.output, output)

    return job
