from typing import TYPE_CHECKING

from cpg_utils import config, hail_batch

if TYPE_CHECKING:
    from hailtop.batch.job import BashJob


def create_filter_wham_jobs(
    input_vcf: str,
    output: str,
) -> 'BashJob':
    in_vcf = hail_batch.get_batch().read_input_group(**{'vcf.gz': input_vcf, 'vcf.gz.tbi': f'{input_vcf}.tbi'})[
        'vcf.gz'
    ]
    job = hail_batch.get_batch().new_bash_job('Filter Wham', attributes={'tool': 'bcftools'})
    job.image(config.config_retrieve(['images', 'bcftools']))
    job.cpu(1).memory('highmem').storage('20Gi')

    job.declare_resource_group(
        output={
            'vcf.bgz': '{root}.vcf.bgz',
            'vcf.bgz.tbi': '{root}.vcf.bgz.tbi',
        }
    )
    job.command(
        'bcftools view -e \'SVTYPE=="DEL" && COUNT(ALGORITHMS)==1 && ALGORITHMS=="wham"\' '
        f'-Oz -o {job.output["vcf.bgz"]} -W=tbi {in_vcf}',
    )
    hail_batch.get_batch().write_output(
        job.output,
        output.replace('.vcf.bgz', ''),
    )
    return job
