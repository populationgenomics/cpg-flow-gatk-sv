"""
All post-batching stages of the GATK-SV workflow
"""

import argparse

from cpg_flow import stage, targets, workflow
from cpg_flow_gatk_sv import utils
from cpg_flow_gatk_sv.jobs import (
    ClusterBatch,
    FilterBatch,
    GatherBatchEvidence,
    GenerateBatchMetrics,
    TrainGCNV,
    MergeBatchSites,
    CombineExclusionLists,
    GenotypeBatch,
    MakeCohortVcf,
)
from cpg_utils import Path, config


@stage.stage
class MakeCohortCombinedPed(stage.CohortStage):
    def expected_outputs(self, cohort: targets.Cohort) -> Path:
        return self.get_stage_cohort_prefix(cohort) / 'ped_with_ref_panel.ped'

    def queue_jobs(self, cohort: targets.Cohort, inputs: stage.StageInput) -> stage.StageOutput:
        output = self.expected_outputs(cohort)

        utils.make_combined_ped(cohort, output)

        return self.make_outputs(target=cohort, data=output)


@stage.stage
class MakeMultiCohortCombinedPed(stage.MultiCohortStage):
    def expected_outputs(self, multicohort: targets.MultiCohort) -> Path:
        return self.prefix / 'ped_with_ref_panel.ped'

    def queue_jobs(self, multicohort: targets.MultiCohort, inputs: stage.StageInput) -> stage.StageOutput:
        # check that there are no overlapping cohorts
        utils.check_for_cohort_overlaps(multicohort)

        output = self.expected_outputs(multicohort)

        utils.make_combined_ped(multicohort, output)

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
                self.get_stage_cohort_prefix(cohort) / f'cohort_gcnv_model_{idx}.tar.gz'
                for idx in range(config.config_retrieve(['workflow', 'model_tar_cnt']))
            ],
        }

    def queue_jobs(self, cohort: targets.Cohort, inputs: stage.StageInput) -> stage.StageOutput:
        outputs = self.expected_outputs(cohort)

        jobs = TrainGCNV.add_train_gcnv_jobs(cohort, output_dict=outputs)

        return self.make_outputs(cohort, data=outputs, jobs=jobs)


@stage.stage(required_stages=[MakeCohortCombinedPed, TrainGcnvStage])
class GatherBatchEvidenceStage(stage.CohortStage):
    """
    https://github.com/broadinstitute/gatk-sv#gather-batch-evidence
    https://github.com/broadinstitute/gatk-sv/blob/master/wdl/GatherBatchEvidence.wdl

    it's critical to separate the ending with a dot, e.g.: `*.sr.txt.gz`,
    These files are passed to `gatk PrintSVEvidence`, that determines file
    format based on the file name.
    It would strongly expect the files to end exactly with either
    `.sr.txt.gz`, `.pe.txt.gz`, or `.sd.txt.gz`, otherwise it would fail with
    "A USER ERROR has occurred: Cannot read file:///cromwell_root/... because
    no suitable codecs found".
    """

    def expected_outputs(self, cohort: targets.Cohort) -> dict[str, Path]:
        """create the output paths for GatherBatchEvidence"""
        ending_by_key = {
            'cnmops_dup': 'DUP.header.bed.gz',
            'cnmops_dup_index': 'DUP.header.bed.gz.tbi',
            'cnmops_del': 'DEL.header.bed.gz',
            'cnmops_del_index': 'DEL.header.bed.gz.tbi',
            'cnmops_large_del': 'DEL.large.bed.gz',
            'cnmops_large_del_index': 'DEL.large.bed.gz.tb',
            'cnmops_large_dup': 'DUP.large.bed.gz',
            'cnmops_large_dup_index': 'DUP.large.bed.gz.tbi',
            'merged_SR': f'{self.name}.sr.txt.gz',
            'merged_SR_index': f'{self.name}.sr.txt.gz.tbi',
            'merged_PE': f'{self.name}.pe.txt.gz',
            'merged_PE_index': f'{self.name}.pe.txt.gz.tbi',
            'merged_BAF': f'{self.name}.baf.txt.gz',
            'merged_BAF_index': f'{self.name}.baf.txt.gz.tbi',
            'merged_bincov': f'{self.name}.RD.txt.gz',
            'merged_bincov_index': f'{self.name}.RD.txt.gz.tbi',
            'median_cov': 'medianCov.transposed.bed',
            'merged_dels': 'DEL.bed.gz',
            'merged_dups': 'DUP.bed.gz',
        }

        # we don't run metrics as standard, only expect the output if we choose to run
        if config.config_retrieve(['resource_overrides', self.name, 'run_matrix_qc'], False):
            ending_by_key.update(
                {
                    'Matrix_QC_plot': '00_matrix_FC_QC.png',
                    'SR_stats': 'SR.QC_matrix.txt',
                    'PE_stats': 'PE.QC_matrix.txt',
                    'BAF_stats': 'BAF.QC_matrix.txt',
                    'RD_stats': 'RD.QC_matrix.txt',
                },
            )

        for caller in utils.SV_CALLERS:
            ending_by_key[f'std_{caller}_vcf_tar'] = f'{caller}.tar.gz'

        return {key: self.get_stage_cohort_prefix(cohort) / fname for key, fname in ending_by_key.items()}

    def queue_jobs(self, cohort: targets.Cohort, inputs: stage.StageInput) -> stage.StageOutput:
        """Add jobs to Batch"""
        train_gcnv_outputs = inputs.as_dict(target=cohort, stage=TrainGcnvStage)

        outputs = self.expected_outputs(cohort)

        jobs = GatherBatchEvidence.submit_gatherbatchevidence_jobs(
            pedigree_input=inputs.as_path(target=cohort, stage=MakeCohortCombinedPed),
            cohort=cohort,
            train_gcnv=train_gcnv_outputs,
            outputs=outputs,
        )

        return self.make_outputs(cohort, data=outputs, jobs=jobs)


@stage.stage(required_stages=[MakeCohortCombinedPed, GatherBatchEvidenceStage])
class ClusterBatchStage(stage.CohortStage):
    """
    https://github.com/broadinstitute/gatk-sv#clusterbatch
    """

    def expected_outputs(self, cohort: targets.Cohort) -> dict:
        """
        * Clustered SV VCFs
        * Clustered depth-only call VCF
        """

        ending_by_key = {}

        for caller in [*utils.SV_CALLERS, 'depth']:
            ending_by_key[f'clustered_{caller}_vcf'] = f'clustered-{caller}.vcf.gz'
            ending_by_key[f'clustered_{caller}_vcf_index'] = f'clustered-{caller}.vcf.gz.tbi'

        return {key: self.get_stage_cohort_prefix(cohort) / fname for key, fname in ending_by_key.items()}

    def queue_jobs(self, cohort: targets.Cohort, inputs: stage.StageInput) -> stage.StageOutput:
        """
        Standardized call VCFs (GatherBatchEvidence)
        Depth-only (DEL/DUP) calls (GatherBatchEvidence)
        """
        batch_evidence_outputs = inputs.as_dict(cohort, GatherBatchEvidenceStage)
        pedigree_input = inputs.as_str(target=cohort, stage=MakeCohortCombinedPed)

        outputs = self.expected_outputs(cohort)

        jobs = ClusterBatch.create_cluster_batch_jobs(
            pedigree_input=pedigree_input,
            cohort=cohort,
            batch_evidence_outputs=batch_evidence_outputs,
            outputs=outputs,
        )
        return self.make_outputs(cohort, data=outputs, jobs=jobs)


@stage.stage(required_stages=[MakeCohortCombinedPed, ClusterBatchStage, GatherBatchEvidenceStage])
class GenerateBatchMetricsStage(stage.CohortStage):
    """
    Generates variant metrics for filtering.
    """

    def expected_outputs(self, cohort: targets.Cohort) -> dict[str, Path]:
        """
        Metrics files
        """

        return {
            'metrics': self.get_stage_cohort_prefix(cohort) / 'metrics.tsv',
            'metrics_common': self.get_stage_cohort_prefix(cohort) / 'metrics_common.tsv',
        }

    def queue_jobs(self, cohort: targets.Cohort, inputs: stage.StageInput) -> stage.StageOutput:
        clusterbatch_outputs = inputs.as_dict(cohort, ClusterBatchStage)
        gatherbatchevidence_outputs = inputs.as_dict(cohort, GatherBatchEvidenceStage)
        pedigree_input = inputs.as_str(target=cohort, stage=MakeCohortCombinedPed)

        outputs = self.expected_outputs(cohort)

        jobs = GenerateBatchMetrics.create_generate_batch_metrics_jobs(
            pedigree_input=pedigree_input,
            cohort=cohort,
            gather_batch_evidence_outputs=gatherbatchevidence_outputs,
            cluster_batch_outputs=clusterbatch_outputs,
            outputs=outputs,
        )

        return self.make_outputs(cohort, data=outputs, jobs=jobs)


@stage.stage(required_stages=[MakeCohortCombinedPed, GenerateBatchMetrics, ClusterBatch])
class FilterBatchStage(stage.CohortStage):
    """
    Filters poor quality variants and filters outlier samples.
    """

    def expected_outputs(self, cohort: targets.Cohort) -> dict:
        """
        * Filtered SV (non-depth-only a.k.a. "PESR") VCF with outlier samples excluded
        * Filtered depth-only call VCF with outlier samples excluded
        * Random forest cutoffs file
        * PED file with outlier samples excluded
        """

        ending_by_key: dict = {
            'filtered_pesr_vcf': 'filtered_pesr_merged.vcf.gz',
            'cutoffs': 'cutoffs',
            'scores': 'updated_scores',
            'RF_intermediate_files': 'RF_intermediate_files.tar.gz',
            'outlier_samples_excluded_file': 'outliers.samples.list',
            'batch_samples_postOutlierExclusion_file': 'outliers_excluded.samples.list',
        }

        for caller in [*utils.SV_CALLERS, 'depth']:
            ending_by_key[f'filtered_{caller}_vcf'] = f'filtered-{caller}.vcf.gz'

            # unsure why, scramble doesn't export this file
            if caller != 'scramble':
                ending_by_key[f'sites_filtered_{caller}_vcf'] = f'sites-filtered-{caller}.vcf.gz'

        ending_by_key['sv_counts'] = [f'{caller}.with_evidence.svcounts.txt' for caller in [*utils.SV_CALLERS, 'depth']]
        ending_by_key['sv_count_plots'] = [
            f'{caller}.with_evidence.all_SVTYPEs.counts_per_sample.png' for caller in [*utils.SV_CALLERS, 'depth']
        ]

        d = {}
        for key, ending in ending_by_key.items():
            if isinstance(ending, str):
                d[key] = self.get_stage_cohort_prefix(cohort) / ending
            elif isinstance(ending, list):
                d[key] = [self.get_stage_cohort_prefix(cohort) / e for e in ending]
        return d

    def queue_jobs(self, cohort: targets.Cohort, inputs: stage.StageInput) -> stage.StageOutput:
        metrics_d = inputs.as_dict(cohort, GenerateBatchMetricsStage)
        clusterbatch_d = inputs.as_dict(cohort, ClusterBatchStage)
        pedigree_input = inputs.as_str(target=cohort, stage=MakeCohortCombinedPed)

        outputs = self.expected_outputs(cohort)

        jobs = FilterBatch.create_filterbatch_jobs(
            pedigree_input=pedigree_input,
            cohort=cohort,
            batch_metrics_outputs=metrics_d,
            cluster_batch_outputs=clusterbatch_d,
            outputs=outputs,
        )

        return self.make_outputs(cohort, data=outputs, jobs=jobs)


@stage.stage(required_stages=FilterBatchStage)
class MergeBatchSitesStage(stage.MultiCohortStage):
    """
    This Stage runs between GenerateBatchMetrics and FilterBatch
    This takes the component VCFs from individual batches and merges them into
    a single VCF (one for PESR and one for depth-only calls).å
    """

    def expected_outputs(self, multicohort: targets.MultiCohort) -> dict[str, Path]:
        return {
            'cohort_pesr_vcf': self.prefix / 'cohort_pesr.vcf.gz',
            'cohort_depth_vcf': self.prefix / 'cohort_depth.vcf.gz',
        }

    def queue_jobs(self, multicohort: targets.MultiCohort, inputs: stage.StageInput) -> stage.StageOutput:
        # take from previous per-cohort outputs
        filter_batch_outputs = inputs.as_dict_by_target(FilterBatchStage)

        pesr_vcfs = [filter_batch_outputs[cohort.id]['filtered_pesr_vcf'] for cohort in multicohort.get_cohorts()]
        depth_vcfs = [filter_batch_outputs[cohort.id]['filtered_depth_vcf'] for cohort in multicohort.get_cohorts()]

        outputs = self.expected_outputs(multicohort)

        jobs = MergeBatchSites.create_mergebatchsites_jobs(
            multicohort=multicohort,
            pesr_vcfs=pesr_vcfs,
            depth_vcfs=depth_vcfs,
            outputs=outputs,
        )

        return self.make_outputs(multicohort, data=outputs, jobs=jobs)


@stage.stage(analysis_type='sv', required_stages=FilterBatchStage)
class CombineExclusionListsStage(stage.MultiCohortStage):
    """
    Takes the per-batch lists of excluded sample IDs and combines
    them into a single file for use in the SV pipeline

    This will be used to remove any filtered samples from consideration in
    subsequent stages, and to remove the CPG ID registration
    """

    def expected_outputs(self, multicohort: targets.MultiCohort) -> Path:
        """
        Create dictionary of names -> output paths
        This variable is a Path to make sure it gets existence checked
        We need this quick stage to run each time
        """

        return self.prefix / 'combined_exclusion_list.txt'

    def queue_jobs(self, multicohort: targets.MultiCohort, inputs: stage.StageInput) -> stage.StageOutput:
        """
        queue job to combine exclusion lists
        """

        filter_batch_outputs = inputs.as_dict_by_target(FilterBatchStage)
        all_filter_lists = [
            str(filter_batch_outputs[cohort.id]['outlier_samples_excluded_file'])
            for cohort in multicohort.get_cohorts()
        ]

        output = self.expected_outputs(multicohort)

        job = CombineExclusionLists.create_combine_exclusion_lists_job(all_filter_lists, output)

        return self.make_outputs(multicohort, data=output, jobs=job)


@stage.stage(
    required_stages=[FilterBatchStage, GatherBatchEvidenceStage, MergeBatchSitesStage],
    analysis_type='sv',
    analysis_keys=['genotyped_depth_vcf', 'genotyped_pesr_vcf'],
)
class GenotypeBatchStage(stage.CohortStage):
    """
    The final Cohort Stage - executed for each individual Batch
    Genotypes a batch of samples across filtered variants combined across all batches.
    Include the additional config file:
    - configs/gatk_sv/use_for_all_workflows.toml; contains all required images and references
    The final CohortStage in this workflow
    """

    def expected_outputs(self, cohort: targets.Cohort) -> dict[str, Path]:
        """
        Filtered SV (non-depth-only a.k.a. "PESR") VCF with outlier samples excluded
        Filtered depth-only call VCF with outlier samples excluded
        PED file with outlier samples excluded
        List of SR pass variants
        List of SR fail variants
        """

        ending_by_key = {
            'sr_bothside_pass': 'genotype_SR_part2_bothside_pass.txt',
            'sr_background_fail': 'genotype_SR_part2_background_fail.txt',
            'trained_PE_metrics': 'pe_metric_file.txt',
            'trained_SR_metrics': 'sr_metric_file.txt',
            'regeno_coverage_medians': 'regeno.coverage_medians_merged.bed',
        }

        for mode in ['pesr', 'depth']:
            ending_by_key |= {
                f'trained_genotype_{mode}_pesr_sepcutoff': f'{mode}.pesr_sepcutoff.txt',
                f'trained_genotype_{mode}_depth_sepcutoff': f'{mode}.depth_sepcutoff.txt',
                f'genotyped_{mode}_vcf': f'{mode}.vcf.gz',
                f'genotyped_{mode}_vcf_index': f'{mode}.vcf.gz.tbi',
            }

        output_hash = workflow.get_workflow().output_version

        return {key: self.get_stage_cohort_prefix(cohort) / output_hash / fname for key, fname in ending_by_key.items()}

    def queue_jobs(self, cohort: targets.Cohort, inputs: stage.StageInput) -> stage.StageOutput:
        filterbatch_outputs = inputs.as_dict(cohort, FilterBatchStage)
        batchevidence_outputs = inputs.as_dict(cohort, GatherBatchEvidenceStage)
        mergebatch_outputs = inputs.as_dict(workflow.get_multicohort(), MergeBatchSitesStage)

        outputs = self.expected_outputs(cohort)

        jobs = GenotypeBatch.create_genotypebatch_jobs(
            cohort=cohort,
            gatherbatchevidence_outputs=batchevidence_outputs,
            filterbatch_outputs=filterbatch_outputs,
            mergebatchsites_outputs=mergebatch_outputs,
            outputs=outputs,
        )

        return self.make_outputs(cohort, data=outputs, jobs=jobs)


@stage.stage(
    required_stages=[
        MakeMultiCohortCombinedPed,
        GatherBatchEvidenceStage,
        GenotypeBatchStage,
        FilterBatchStage,
    ]
)
class MakeCohortVcfStage(stage.MultiCohortStage):
    """
    Combines variants across multiple batches, resolves complex variants, re-genotypes, and performs final VCF clean-up.
    """

    def expected_outputs(self, multicohort: targets.MultiCohort) -> dict:
        """create output paths"""
        outputs = {
            'vcf': self.prefix / 'cleaned.vcf.gz',
            'vcf_index': self.prefix / 'cleaned.vcf.gz.tbi',
        }
        if config.config_retrieve(['MakeCohortVcf', 'MainVcfQc', 'do_per_sample_qc'], False):
            outputs |= {'vcf_qc': self.prefix / 'cleaned_SV_VCF_QC_output.tar.gz'}

        return outputs

    def queue_jobs(self, multicohort: stage.MultiCohort, inputs: stage.StageInput) -> stage.StageOutput:
        """
        Instead of taking a direct dependency on the previous stage,
        we use the output hash to find all the previous batches
        """
        gatherbatchevidence_outputs = inputs.as_dict_by_target(GatherBatchEvidenceStage)
        genotypebatch_outputs = inputs.as_dict_by_target(GenotypeBatchStage)
        filterbatch_outputs = inputs.as_dict_by_target(FilterBatchStage)
        pedigree_input = inputs.as_path(target=multicohort, stage=MakeMultiCohortCombinedPed)

        outputs = self.expected_outputs(multicohort)

        jobs = MakeCohortVcf.create_makecohortvcf_jobs(
            multicohort,
            pedigree_input,
            gatherbatchevidence_outputs,
            genotypebatch_outputs,
            filterbatch_outputs,
            outputs,
        )

        return self.make_outputs(multicohort, data=outputs, jobs=jobs)


def cli_main():
    """
    CLI entrypoint - starts up the workflow
    """
    parser = argparse.ArgumentParser()
    parser.add_argument('--dry_run', action='store_true', help='Dry run')
    args = parser.parse_args()

    workflow.run_workflow(
        stages=[
            TrainGcnvStage,
            MakeCohortCombinedPed,
            MakeMultiCohortCombinedPed,
            GatherBatchEvidenceStage,
            ClusterBatch,
            GenerateBatchMetricsStage,
            FilterBatchStage,
            MergeBatchSitesStage,
            CombineExclusionListsStage,
            GenotypeBatchStage,
            MakeCohortVcfStage,
        ],
        dry_run=args.dry_run,
    )


if __name__ == '__main__':
    cli_main()
