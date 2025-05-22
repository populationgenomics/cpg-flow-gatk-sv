""" """

from typing import TYPE_CHECKING

from cpg_flow import targets
from cpg_utils import Path
from cpg_utils.config import AR_GUID_NAME, try_get_ar_guid

from cpg_flow_gatk_sv import utils

if TYPE_CHECKING:
    from hailtop.batch.job import BashJob


def create_evidence_qc_jobs(
    input_dict: dict[str, Path],
    output_dict: dict[str, Path],
    cohort: targets.Cohort,
) -> list['BashJob']:
    """

    Args:
        input_dict ():
        output_dict ():
        cohort ():

    Returns:

    """

    sgids = cohort.get_sequencing_group_ids()

    cromwell_input_dict: dict = {
        'batch': cohort.id,
        'samples': sgids,
        'run_vcf_qc': True,
        'counts': [str(input_dict[sid]['coverage_counts']) for sid in sgids],
    }

    for caller in utils.SV_CALLERS:
        cromwell_input_dict[f'{caller}_vcfs'] = [str(input_dict[sid][f'{caller}_vcf']) for sid in sgids]

    cromwell_input_dict |= utils.get_images(
        ['sv_base_mini_docker', 'sv_base_docker', 'sv_pipeline_docker', 'sv_pipeline_qc_docker'],
    )

    cromwell_input_dict |= utils.get_references(['genome_file', 'wgd_scoring_mask'])

    # runs for approx 5 hours, depending on sample count
    return utils.add_gatk_sv_jobs(
        dataset=cohort.dataset,
        wfl_name='EvidenceQC',
        input_dict=cromwell_input_dict,
        expected_out_dict=output_dict,
        labels={'stage': 'evidenceqc', AR_GUID_NAME: try_get_ar_guid()},
        job_size=utils.CromwellJobSizes.MEDIUM,
    )
