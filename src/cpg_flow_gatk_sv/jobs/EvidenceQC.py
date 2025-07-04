""" """

from typing import TYPE_CHECKING

from cpg_flow import targets
from cpg_flow_gatk_sv import utils
from cpg_utils import Path, config

if TYPE_CHECKING:
    from hailtop.batch.job import BashJob


def create_evidence_qc_jobs(
    input_dict: dict[str, Path],
    output_dict: dict[str, Path],
    cohort: targets.Cohort,
) -> list['BashJob']:
    sgids = cohort.get_sequencing_group_ids()

    cromwell_input_dict: dict = {
        'batch': cohort.id,
        'samples': sgids,
        'run_vcf_qc': True,
        'counts': [str(input_dict[sid]['coverage_counts']) for sid in sgids],
        'sv_base_docker': config.config_retrieve(['images', 'sv_base_docker']),
        'sv_pipeline_qc_docker': config.config_retrieve(['images', 'sv_pipeline_qc_docker']),
        'sv_base_mini_docker': config.config_retrieve(['images', 'sv_base_mini_docker']),
        'sv_pipeline_docker': config.config_retrieve(['images', 'sv_pipeline_docker']),
        'genome_file': config.config_retrieve(['references', 'wgd_scoring_mask']),
    }

    for caller in utils.SV_CALLERS:
        cromwell_input_dict[f'{caller}_vcfs'] = [str(input_dict[sid][f'{caller}_vcf']) for sid in sgids]

    # runs for approx 5 hours, depending on sample count
    return utils.add_gatk_sv_jobs(
        dataset=cohort.dataset,
        wfl_name='EvidenceQC',
        input_dict=cromwell_input_dict,
        expected_out_dict=output_dict,
        labels={'stage': 'evidenceqc', config.AR_GUID_NAME: config.try_get_ar_guid()},
        job_size=utils.CromwellJobSizes.MEDIUM,
    )
