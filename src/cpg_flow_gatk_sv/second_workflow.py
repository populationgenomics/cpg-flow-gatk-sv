"""
All post-batching stages of the GATK-SV workflow
"""

from cpg_utils import Path
from cpg_flow import stage, targets, workflow
from cpg_flow_gatk_sv.utils import check_for_cohort_overlaps, make_combined_ped


@stage.stage
class MakeCohortCombinedPed(stage.CohortStage):
    def expected_outputs(self, cohort: targets.Cohort) -> Path:
        return self.get_stage_cohort_prefix(cohort) / 'ped_with_ref_panel.ped'

    def queue_jobs(self, cohort: targets.Cohort, inputs: stage.StageInput) -> stage.StageOutput:
        output = self.expected_outputs(cohort)

        make_combined_ped(cohort, output)

        return self.make_outputs(target=cohort, data=output)


@stage.stage
class MakeMultiCohortCombinedPed(stage.MultiCohortStage):
    def expected_outputs(self, multicohort: targets.MultiCohort) -> Path:
        return self.prefix / 'ped_with_ref_panel.ped'

    def queue_jobs(self, multicohort: targets.MultiCohort, inputs: stage.StageInput) -> stage.StageOutput:
        # check that there are no overlapping cohorts
        check_for_cohort_overlaps(multicohort)

        output = self.expected_outputs(multicohort)

        make_combined_ped(multicohort, output)

        return self.make_outputs(target=multicohort, data=output)
