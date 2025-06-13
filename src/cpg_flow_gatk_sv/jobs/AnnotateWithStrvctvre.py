from typing import TYPE_CHECKING

from cpg_flow_gatk_sv import utils
from cpg_utils import config, hail_batch

if TYPE_CHECKING:
    from hailtop.batch.job import BashJob


def create_strvctvre_jobs(
    input_vcf: str,
    output: str,
    name: str = 'AnnotateVcfWithStrvctvre',
) -> 'BashJob':
    input_vcf = hail_batch.get_batch().read_input_group(
        vcf=input_vcf,
        vcf_index=f'{input_vcf}.tbi',
    )['vcf']

    job = hail_batch.get_batch().new_bash_job('StrVCTVRE', {'tool': 'strvctvre'})

    job.image(config.config_retrieve(['images', 'strvctvre']))
    job.storage(config.config_retrieve(['resource_overrides', name, 'storage'], '10Gi'))
    job.memory(config.config_retrieve(['resource_overrides', name, 'memory'], '16Gi'))

    strvctvre_phylop = utils.get_references(['strvctvre_phylop'])['strvctvre_phylop']

    local_phylop = hail_batch.get_batch().read_input(strvctvre_phylop)

    job.declare_resource_group(
        output={
            'vcf.gz': '{root}.vcf.gz',
            'vcf.gz.tbi': '{root}.vcf.gz.tbi',
        }
    )

    # run strvctvre
    job.command(
        f'python StrVCTVRE.py -i {input_vcf} -o {job.output["vcf.gz"]} -f vcf -p {local_phylop}',
    )
    job.command(f'tabix {job.output["vcf.gz"]}')

    hail_batch.get_batch().write_output(job.output, output.replace('.vcf.gz', ''))
    return job
