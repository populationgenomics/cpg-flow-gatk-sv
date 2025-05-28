from typing import TYPE_CHECKING

from cpg_utils import config, hail_batch

if TYPE_CHECKING:
    from hailtop.batch.job import BashJob


def create_spicy_jobs(
    input_vcf: str,
    skip_prior_names: bool,
    output: str,
) -> 'BashJob':
    local_input_vcf = hail_batch.get_batch().read_input(input_vcf)

    # update the IDs using a PythonJob
    job = hail_batch.get_batch().new_bash_job('rename_sv_ids')
    job.storage('10Gi')
    job.image(config.config_retrieve(['workflow', 'driver_image']))

    job.command(f"""
    python -m cpg_flow_gatk_sv.scripts.rename_sv_ids \
        --input_vcf {local_input_vcf} \
        --output_vcf {job.output} \
        {'--skip_prior_names' if skip_prior_names else ''} \
    """)

    # then compress & run tabix on that plain text result
    bcftools_job = hail_batch.get_batch().new_job('bgzip and tabix')
    bcftools_job.image(config.config_retrieve(['images', 'bcftools']))
    bcftools_job.declare_resource_group(
        output={
            'vcf.bgz': '{root}.vcf.bgz',
            'vcf.bgz.tbi': '{root}.vcf.bgz.tbi',
        }
    )
    bcftools_job.command(f'bcftools view {job.output} -Oz -W=tbi -o {bcftools_job.output["vcf.bgz"]}')

    # get the output root to write to
    hail_batch.get_batch().write_output(bcftools_job.output, output.removesuffix('.vcf.bgz'))
    return [job, bcftools_job]
