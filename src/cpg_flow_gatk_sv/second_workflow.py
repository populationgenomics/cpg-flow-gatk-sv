"""
All post-batching stages of the GATK-SV workflow
"""

import argparse

from cpg_utils import Path
from cpg_flow import stage, targets, workflow
from cpg_flow_gatk_sv.utils import check_for_cohort_overlaps, make_combined_ped
from cpg_flow_gatk_sv.jobs import TrainGCNV


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


@stage.stage(required_stages=[MakeCohortCombinedPed])
class TrainGcnvStage(stage.CohortStage):
    """
    Runs TrainGCNV to generate a model trained on a selection of samples within this Cohort, used in GatherBatchEvidence
    run https://github.com/broadinstitute/gatk-sv/blob/main/wdl/TrainGCNV.wdl
    config: https://github.com/broadinstitute/gatk-sv/blob/main/inputs/templates/terra_workspaces/cohort_mode/workflow_configurations/TrainGCNV.json.tmpl
    """

    def expected_outputs(self, cohort: targets.Cohort) -> dict[str, Path]:
        return {
            'cohort_contig_ploidy_model_tar': self.get_stage_cohort_prefix(cohort)
            / 'cohort_contig_ploidy_model.tar.gz',
            'cohort_contig_ploidy_calls_tar': self.get_stage_cohort_prefix(cohort)
            / 'cohort_contig_ploidy_calls.tar.gz',
            'cohort_gcnv_model_tars': [
                self.get_stage_cohort_prefix(cohort) / f'cohort_gcnv_model_{idx}.tar.gz' for idx in range(290)
            ],  # I have absolutely no idea how long this will be
        }

    def queue_jobs(self, cohort: targets.Cohort, inputs: stage.StageInput) -> stage.StageOutput:
        outputs = self.expected_outputs(cohort)

        jobs = TrainGCNV.add_train_gcnv_jobs(cohort, output_dict=outputs)

        return self.make_outputs(cohort, data=outputs, jobs=jobs)


def cli_main():
    """
    CLI entrypoint - starts up the workflow
    """
    parser = argparse.ArgumentParser()
    parser.add_argument('--dry_run', action='store_true', help='Dry run')
    args = parser.parse_args()

    workflow.run_workflow(
        stages=[TrainGcnvStage, MakeCohortCombinedPed, MakeMultiCohortCombinedPed],
        dry_run=args.dry_run,
    )


if __name__ == '__main__':
    cli_main()
