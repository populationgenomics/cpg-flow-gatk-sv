import os
from typing import TYPE_CHECKING, Any

from cpg_flow import targets
from cpg_flow_gatk_sv import utils
from cpg_utils import Path, config, hail_batch, to_path

if TYPE_CHECKING:
    from hailtop.batch.job import BashJob


STAGE_NAME: str = 'GatherSampleEvidence'


def create_gather_sample_evidence_jobs(
    sg: targets.SequencingGroup,
    expected_outputs: dict[str, Path],
) -> list['BashJob']:
    fasta_file = config.config_retrieve(['workflow', 'ref_fasta'])

    input_dict: dict[str, Any] = {
        'bam_or_cram_file': str(sg.cram),
        'bam_or_cram_index': str(sg.cram) + '.crai',
        'cloud_sdk_docker': config.config_retrieve(['images', 'cloud_sdk_docker']),
        'gatk_docker': config.config_retrieve(['images', 'gatk_docker']),
        'gatk_docker_pesr_override': config.config_retrieve(['images', 'gatk_docker_pesr_override']),
        'genomes_in_the_cloud_docker': config.config_retrieve(['images', 'genomes_in_the_cloud_docker']),
        'manta_docker': config.config_retrieve(['images', 'manta_docker']),
        'manta_region_bed': config.config_retrieve(['references', 'manta_region_bed']),
        'manta_region_bed_index': config.config_retrieve(['references', 'manta_region_bed_index']),
        'mei_bed': config.config_retrieve(['references', 'mei_bed']),
        # a cost-improvement in cloud environments
        'move_bam_or_cram_files': True,
        'preprocessed_intervals': config.config_retrieve(['references', 'preprocessed_intervals']),
        'primary_contigs_fai': config.config_retrieve(['references', 'primary_contigs_fai']),
        'primary_contigs_list': config.config_retrieve(['references', 'primary_contigs_list']),
        'reference_dict': to_path(fasta_file).with_suffix('.dict'),
        'reference_fasta': fasta_file,
        'reference_index': f'{fasta_file}.fai',
        'reference_version': '38',
        'sample_id': sg.id,
        'samtools_cloud_docker': config.config_retrieve(['images', 'samtools_cloud_docker']),
        'scramble_docker': config.config_retrieve(['images', 'scramble_docker']),
        'sd_locs_vcf': config.config_retrieve(['references', 'dbsnp_vcf']),
        'sv_base_mini_docker': config.config_retrieve(['images', 'sv_base_mini_docker']),
        'sv_pipeline_docker': config.config_retrieve(['images', 'sv_pipeline_docker']),
        'wham_docker': config.config_retrieve(['images', 'wham_docker']),
        'wham_include_list_bed_file': config.config_retrieve(['references', 'wham_include_list_bed_file']),
    }

    # If DRAGEN input is going to be used, first the input parameter 'is_dragen_3_7_8' needs to be set to True
    # then some parameters need to be added to the input_dict to enable BWA to be run

    if only_jobs := config.config_retrieve(['workflow', STAGE_NAME, 'only_jobs'], None):
        # if only_jobs is set, only run the specified jobs
        # this is useful for samples which need to re-run specific jobs
        # e.g. if manta failed and needs to be re-run with more memory

        # disable the evidence collection jobs if they're not in only_jobs
        if 'coverage_counts' not in only_jobs:
            input_dict['collect_coverage'] = False
        if 'pesr' not in only_jobs:
            input_dict['collect_pesr'] = False

        if 'scramble' in only_jobs and 'manta' not in only_jobs:
            # if Scramble is being run, but Manta is not, manta_vcf and index becomes a required input
            if to_path(expected_outputs['manta_vcf']).exists():
                # if the Manta VCF exists, use it as input and remove it from expected outputs
                input_dict['manta_vcf_input'] = expected_outputs.pop('manta_vcf')
                input_dict['manta_vcf_index_input'] = expected_outputs.pop('manta_index')
            else:
                # if the Manta VCF does not exist, run Manta as well
                only_jobs.append('manta')

        # disable the caller jobs that are not in only_jobs by nulling their docker image
        for key, val in input_dict.items():
            if key in [f'{caller}_docker' for caller in utils.SV_CALLERS]:
                caller = key.removesuffix('_docker')
                input_dict[key] = val if caller in only_jobs else None

    # billing labels!
    # https://cromwell.readthedocs.io/en/stable/wf_options/Google/
    # these must conform to the regex [a-z]([-a-z0-9]*[a-z0-9])?
    billing_labels = {
        'dataset': sg.dataset.name,  # already lowercase
        'sequencing-group': sg.id.lower(),
        'stage': STAGE_NAME.lower(),
        config.AR_GUID_NAME: config.try_get_ar_guid(),
    }

    # add some max-polling interval jitter for each sample
    # cromwell_status_poll_interval is a number (int, seconds)
    # this is used to determine how often to poll Cromwell for completion status
    # we alter the per-sample maximum to be between 5 and 30 minutes for this
    # long-running job, so samples poll on different intervals, spreading load
    gather_sample_evidence_jobs = utils.add_gatk_sv_jobs(
        dataset=sg.dataset,
        wfl_name=STAGE_NAME,
        input_dict=input_dict,
        expected_out_dict=expected_outputs,
        sequencing_group_id=sg.id,
        labels=billing_labels,
        job_size=utils.CromwellJobSizes.LARGE,
    )

    deletion_job = hail_batch.get_batch().new_bash_job(f'Delete GatherSampleEvidence tmp for {sg.id}')
    deletion_job.image(config.config_retrieve(['workflow', 'driver_image']))
    # set explicit dependency - no variant calling success, no deletion
    deletion_job.depends_on(*gather_sample_evidence_jobs)
    tmp_bucket = config.config_retrieve(['storage', 'default', 'tmp'])
    delete_path = os.path.join(tmp_bucket, 'cromwell', STAGE_NAME, sg.id)
    deletion_job.command(f'gcloud storage rm -r {delete_path}')

    return [*gather_sample_evidence_jobs, deletion_job]
