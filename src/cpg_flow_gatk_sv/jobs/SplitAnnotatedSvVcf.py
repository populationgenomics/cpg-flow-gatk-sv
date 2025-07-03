from typing import TYPE_CHECKING

from cpg_flow import targets
from cpg_utils import config, hail_batch

if TYPE_CHECKING:
    from hailtop.batch.job import BashJob


def create_split_vcf_by_dataset_job(
    dataset: targets.Dataset,
    input_vcf: str,
    dataset_sgid_file: str,
    output: str,
    job_attrs: dict[str, str],
) -> 'BashJob':
    """Split the MultiCohort VCF by dataset."""

    local_vcf = hail_batch.get_batch().read_input(input_vcf)
    local_sgids = hail_batch.get_batch().read_input(dataset_sgid_file)

    job = hail_batch.get_batch().new_bash_job(
        name=f'SplitAnnotatedSvVcfByDataset: {dataset.name}',
        attributes=job_attrs | {'tool': 'bcftools'},
    )
    job.image(config.config_retrieve(['images', 'bcftools']))
    job.cpu(1).memory('highmem').storage('10Gi')
    job.declare_resource_group(
        output={
            'vcf.bgz': '{root}.vcf.bgz',
            'vcf.bgz.tbi': '{root}.vcf.bgz.tbi',
        }
    )

    job.command(
        f"""
        bcftools view \\
            --force-samples \\
            -S {local_sgids} \\
            -Oz \\
            -o {job.output['vcf.bgz']} \\
            --W=tbi \\
            {local_vcf}
        """,
    )

    hail_batch.get_batch().write_output(job.output, output.removesuffix('.vcf.bgz'))

    return job
