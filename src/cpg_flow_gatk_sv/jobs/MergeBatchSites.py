from typing import TYPE_CHECKING

from cpg_flow import targets
from cpg_flow_gatk_sv import utils
from cpg_utils import Path, config

if TYPE_CHECKING:
    from hailtop.batch.job import BashJob


def create_mergebatchsites_jobs(
    multicohort: targets.MultiCohort,
    pesr_vcfs: list[Path],
    depth_vcfs: list[Path],
    outputs: dict[str, Path],
) -> list['BashJob']:
    input_dict = {
        'cohort': multicohort.name,
        'depth_vcfs': depth_vcfs,
        'pesr_vcfs': pesr_vcfs,
        'sv_pipeline_docker': config.config_retrieve(['images', 'sv_pipeline_docker']),
    }

    return utils.add_gatk_sv_jobs(
        dataset=multicohort.analysis_dataset,
        wfl_name='MergeBatchSites',
        input_dict=input_dict,
        expected_out_dict=outputs,
        labels={'stage': 'mergebatchsites', config.AR_GUID_NAME: config.try_get_ar_guid()},
    )
