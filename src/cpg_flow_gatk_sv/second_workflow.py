"""
All post-batching stages of the GATK-SV workflow
"""

from cpg_utils import Path
from cpg_flow import stage, targets, workflow
from cpg_flow_gatk_sv.utils import make_combined_ped


@stage.stage
class MakeCohortCombinedPed(stage.CohortStage):
    def expected_outputs(self, cohort: targets.Cohort) -> Path:
        return self.get_stage_cohort_prefix(cohort) / 'ped_with_ref_panel.ped'

    def queue_jobs(self, cohort: targets.Cohort, inputs: stage.StageInput) -> stage.StageOutput:
        output = self.expected_outputs(cohort)

        make_combined_ped(cohort, output)

        return self.make_outputs(target=cohort, data=output)
