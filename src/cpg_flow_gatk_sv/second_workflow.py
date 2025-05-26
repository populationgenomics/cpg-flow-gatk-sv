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
    FormatVcfForGatk,
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


@stage.stage(required_stages=[MakeMultiCohortCombinedPed, MakeCohortVcfStage])
class FormatVcfForGatkStage(stage.MultiCohortStage):
    def expected_outputs(self, multicohort: targets.MultiCohort) -> dict:
        return {
            'gatk_formatted_vcf': self.prefix / 'gatk_formatted.vcf.gz',
            'gatk_formatted_vcf_index': self.prefix / 'gatk_formatted.vcf.gz.tbi',
        }

    def queue_jobs(self, multicohort: targets.MultiCohort, inputs: stage.StageInput) -> stage.StageOutput:
        pedigree_input = inputs.as_str(target=multicohort, stage=MakeMultiCohortCombinedPed)
        vcf_input = inputs.as_str(multicohort, MakeCohortVcfStage, 'vcf')
        outputs = self.expected_outputs(multicohort)

        jobs = FormatVcfForGatk.create_formatvcf_jobs(
            multicohort,
            pedigree_input,
            vcf_input,
            outputs,
        )

        return self.make_outputs(multicohort, data=outputs, jobs=jobs)


# @stage(required_stages=[MakeMultiCohortCombinedPed, ClusterBatch, MakeCohortVcf])
# class JoinRawCalls(MultiCohortStage):
#     """
#     Joins all individually clustered caller results
#     """
#
#     def expected_outputs(self, multicohort: MultiCohort) -> dict[str, Path]:
#         return {
#             'joined_raw_calls_vcf': self.prefix / 'raw_clustered_calls.vcf.gz',
#             'joined_raw_calls_vcf_index': self.prefix / 'raw_clustered_calls.vcf.gz.tbi',
#         }
#
#     def queue_jobs(self, multicohort: MultiCohort, inputs: StageInput) -> StageOutput:
#         pedigree_input = inputs.as_path(target=multicohort, stage=MakeMultiCohortCombinedPed, key='multicohort_ped')
#         input_dict: dict[str, Any] = {
#             'FormatVcfForGatk.formatter_args': '--fix-end',
#             'prefix': multicohort.name,
#             'ped_file': str(pedigree_input),
#             'reference_fasta': get_fasta(),
#             'reference_fasta_fai': str(get_fasta()) + '.fai',
#             'reference_dict': str(get_fasta().with_suffix('.dict')),
#         }
#         input_dict |= get_images(['gatk_docker', 'sv_pipeline_docker', 'sv_base_mini_docker'])
#         input_dict |= get_references([{'contig_list': 'primary_contigs_list'}])
#
#         # add all clustered _caller_ files, plus indices
#         clusterbatch_outputs = inputs.as_dict_by_target(ClusterBatch)
#
#         # get the names of all contained cohorts
#         all_batch_names: list[str] = [cohort.id for cohort in multicohort.get_cohorts()]
#         for caller in SV_CALLERS + ['depth']:
#             input_dict[f'clustered_{caller}_vcfs'] = [
#                 clusterbatch_outputs[cohort][f'clustered_{caller}_vcf'] for cohort in all_batch_names
#             ]
#             input_dict[f'clustered_{caller}_vcf_indexes'] = [
#                 clusterbatch_outputs[cohort][f'clustered_{caller}_vcf_index'] for cohort in all_batch_names
#             ]
#
#         expected_d = self.expected_outputs(multicohort)
#
#         jobs = add_gatk_sv_jobs(
#             dataset=multicohort.analysis_dataset,
#             wfl_name=self.name,
#             input_dict=input_dict,
#             expected_out_dict=expected_d,
#             labels={'stage': self.name.lower(), AR_GUID_NAME: try_get_ar_guid()},
#         )
#         return self.make_outputs(multicohort, data=expected_d, jobs=jobs)
#
#
# @stage(required_stages=[JoinRawCalls, FormatVcfForGatk])
# class SVConcordance(MultiCohortStage):
#     """
#     Takes the clean VCF and reformat for GATK intake
#     """
#
#     def expected_outputs(self, multicohort: MultiCohort) -> dict:
#         """
#         create dictionary of names -> output paths
#         """
#
#         return {
#             'concordance_vcf': self.prefix / 'sv_concordance.vcf.gz',
#             'concordance_vcf_index': self.prefix / 'sv_concordance.vcf.gz.tbi',
#         }
#
#     def queue_jobs(self, multicohort: MultiCohort, inputs: StageInput) -> StageOutput | None:
#         """
#         configure and queue jobs for SV concordance
#         """
#
#         input_dict: dict[str, Any] = {
#             'output_prefix': multicohort.name,
#             'reference_dict': str(get_fasta().with_suffix('.dict')),
#             'eval_vcf': inputs.as_dict(multicohort, FormatVcfForGatk)['gatk_formatted_vcf'],
#             'truth_vcf': inputs.as_dict(multicohort, JoinRawCalls)['joined_raw_calls_vcf'],
#         }
#         input_dict |= get_images(['gatk_docker', 'sv_base_mini_docker'])
#         input_dict |= get_references([{'contig_list': 'primary_contigs_list'}])
#
#         expected_d = self.expected_outputs(multicohort)
#
#         jobs = add_gatk_sv_jobs(
#             dataset=multicohort.analysis_dataset,
#             wfl_name=self.name,
#             input_dict=input_dict,
#             expected_out_dict=expected_d,
#             labels={'stage': self.name.lower(), AR_GUID_NAME: try_get_ar_guid()},
#         )
#         return self.make_outputs(multicohort, data=expected_d, jobs=jobs)
#
#
# @stage(required_stages=[MakeMultiCohortCombinedPed, SVConcordance])
# class GeneratePloidyTable(MultiCohortStage):
#     """
#     Quick PythonJob to generate a ploidy table
#     Calls a homebrewed version of this table generator:
#     github.com/broadinstitute/gatk-sv/blob/main/src/sv-pipeline/scripts/ploidy_table_from_ped.py
#     """
#
#     def expected_outputs(self, multicohort: MultiCohort) -> dict[str, Path]:
#         """
#         only one output, the ploidy table
#         """
#
#         return {'ploidy_table': self.prefix / 'ploidy_table.txt'}
#
#     def queue_jobs(self, multicohort: MultiCohort, inputs: StageInput) -> StageOutput | None:
#         pedigree_input = inputs.as_path(target=multicohort, stage=MakeMultiCohortCombinedPed, key='multicohort_ped')
#
#         py_job = get_batch().new_python_job('create_ploidy_table')
#         py_job.image(config_retrieve(['workflow', 'driver_image']))
#
#         contig_path = get_references(['primary_contigs_list'])['primary_contigs_list']
#
#         expected_d = self.expected_outputs(multicohort)
#         py_job.call(
#             ploidy_table_from_ped.generate_ploidy_table,
#             str(pedigree_input),
#             contig_path,
#             str(expected_d['ploidy_table']),
#         )
#
#         return self.make_outputs(multicohort, data=expected_d, jobs=py_job)


# @stage(required_stages=[MakeMultiCohortCombinedPed, GeneratePloidyTable, SVConcordance])
# class FilterGenotypes(MultiCohortStage):
#     """
#     Steps required to post-filter called genotypes
#     """
#
#     def expected_outputs(self, multicohort: MultiCohort) -> dict:
#         """
#         create dictionary of names -> output paths
#         """
#
#         return {
#             'filtered_vcf': self.prefix / 'filtered.vcf.gz',
#             'filtered_vcf_index': self.prefix / 'filtered.vcf.gz.tbi',
#             'unfiltered_recalibrated_vcf': self.prefix / 'unfiltered_recalibrated.vcf.gz',
#             'unfiltered_recalibrated_vcf_index': self.prefix / 'unfiltered_recalibrated.vcf.gz.tbi',
#         }
#
#     def queue_jobs(self, multicohort: MultiCohort, inputs: StageInput) -> StageOutput | None:
#         pedigree_input = inputs.as_path(target=multicohort, stage=MakeMultiCohortCombinedPed, key='multicohort_ped')
#         input_dict = {
#             'output_prefix': multicohort.name,
#             'vcf': inputs.as_dict(multicohort, SVConcordance)['concordance_vcf'],
#             'ploidy_table': inputs.as_dict(multicohort, GeneratePloidyTable)['ploidy_table'],
#             'ped_file': str(pedigree_input),
#             'fmax_beta': config_retrieve(['references', 'gatk_sv', 'fmax_beta'], 0.4),
#             'recalibrate_gq_args': config_retrieve(['references', 'gatk_sv', 'recalibrate_gq_args']),
#             'sl_filter_args': config_retrieve(['references', 'gatk_sv', 'sl_filter_args']),
#         }
#         assert input_dict['recalibrate_gq_args'] and input_dict['sl_filter_args'], input_dict
#
#         input_dict |= get_images(['linux_docker', 'sv_base_mini_docker', 'sv_pipeline_docker'])
#         # use a non-standard GATK image containing required filtering tool
#         input_dict['gatk_docker'] = get_images(['gq_recalibrator_docker'])['gq_recalibrator_docker']
#         input_dict |= get_references(
#             [{'gq_recalibrator_model_file': 'aou_filtering_model'}, 'primary_contigs_fai'],
#         )
#
#         # something a little trickier - we need to get various genome tracks
#         input_dict['genome_tracks'] = list(
#             get_references(config_retrieve(['references', 'gatk_sv', 'genome_tracks'], [])).values(),
#         )
#
#         expected_d = self.expected_outputs(multicohort)
#
#         jobs = add_gatk_sv_jobs(
#             dataset=multicohort.analysis_dataset,
#             wfl_name=self.name,
#             input_dict=input_dict,
#             expected_out_dict=expected_d,
#             labels={'stage': self.name.lower(), AR_GUID_NAME: try_get_ar_guid()},
#         )
#         return self.make_outputs(multicohort, data=expected_d, jobs=jobs)
#
#
# @stage(required_stages=FilterGenotypes)
# class UpdateStructuralVariantIDs(MultiCohortStage):
#     """
#     Runs SVConcordance between the results of this callset and the results of a previous callset
#     This causes the Variant IDs of a matching variant to be updated to the previous callset's ID
#     Consistency of Variant ID is crucial to Seqr/AIP identifying the same variant across different callsets
#     By default GATK-SV creates an auto-incrementing ID, rather than one based on variant attributes
#     If a new call is added at a chromosome start for any sample, the ID for all variants on the chromosome
#     will be shifted up, meaning that Seqr labels/any other process linked to the ID will be incorrect
#     """
#
#     def expected_outputs(self, multicohort: MultiCohort) -> dict[str, Path]:
#         return {
#             'concordance_vcf': self.prefix / 'updated_ids.vcf.gz',
#             'concordance_vcf_index': self.prefix / 'updated_ids.vcf.gz.tbi',
#         }
#
#     def queue_jobs(self, multicohort: MultiCohort, inputs: StageInput) -> StageOutput:
#
#         # allow for no prior name/IDs
#         if not (spicy_vcf := query_for_spicy_vcf(multicohort.analysis_dataset.name)):
#             get_logger().info('No previous Spicy VCF found for {cohort.analysis_dataset.name}')
#             return self.make_outputs(multicohort, skipped=True)
#
#         expected_d = self.expected_outputs(multicohort)
#
#         # run concordance between this and the previous VCF
#         input_dict: dict = {
#             'output_prefix': multicohort.name,
#             'reference_dict': str(get_fasta().with_suffix('.dict')),
#             'eval_vcf': inputs.as_path(multicohort, FilterGenotypes, key='filtered_vcf'),
#             'truth_vcf': spicy_vcf,
#         }
#         input_dict |= get_images(['gatk_docker', 'sv_base_mini_docker'])
#         input_dict |= get_references([{'contig_list': 'primary_contigs_list'}])
#
#         jobs = add_gatk_sv_jobs(
#             dataset=multicohort.analysis_dataset,
#             wfl_name='SVConcordance',
#             input_dict=input_dict,
#             expected_out_dict=expected_d,
#             labels={'stage': self.name.lower(), AR_GUID_NAME: try_get_ar_guid()},
#         )
#         return self.make_outputs(multicohort, data=expected_d, jobs=jobs)
#
#
# @stage(
#     required_stages=[FilterGenotypes, UpdateStructuralVariantIDs],
#     analysis_type='sv',
#     analysis_keys=['wham_filtered_vcf'],
# )
# class FilterWham(MultiCohortStage):
#     """
#     Filters the VCF to remove deletions only called by Wham
#     github.com/broadinstitute/gatk-sv/blob/main/wdl/ApplyManualVariantFilter.wdl
#     """
#
#     def expected_outputs(self, multicohort: MultiCohort) -> dict[str, Path]:
#         # index here is implicit
#         return {'wham_filtered_vcf': self.prefix / 'filtered.vcf.bgz'}
#
#     def queue_jobs(self, multicohort: MultiCohort, inputs: StageInput) -> StageOutput:
#         """
#         configure and queue jobs for SV annotation
#         passing the VCF Index has become implicit, which may be a problem for us
#         """
#         # read the concordance-with-prev-batch VCF if appropriate, otherwise use the filtered VCF
#         if query_for_spicy_vcf(multicohort.analysis_dataset.name):
#             get_logger().info(f'Variant IDs were updated for {multicohort.analysis_dataset.name}')
#             input_vcf = inputs.as_dict(multicohort, UpdateStructuralVariantIDs)['concordance_vcf']
#         else:
#             get_logger().info(f'No Spicy VCF was found, default IDs for {multicohort.analysis_dataset.name}')
#             input_vcf = inputs.as_dict(multicohort, FilterGenotypes)['filtered_vcf']
#
#         in_vcf = get_batch().read_input_group(**{'vcf.gz': input_vcf, 'vcf.gz.tbi': f'{input_vcf}.tbi'})['vcf.gz']
#         job = get_batch().new_job('Filter Wham', attributes={'tool': 'bcftools'})
#         job.image(image_path('bcftools')).cpu(1).memory('highmem').storage('20Gi')
#         job.declare_resource_group(output={'vcf.bgz': '{root}.vcf.bgz', 'vcf.bgz.tbi': '{root}.vcf.bgz.tbi'})
#         job.command(
#             'bcftools view -e \'SVTYPE=="DEL" && COUNT(ALGORITHMS)==1 && ALGORITHMS=="wham"\' '
#             f'{in_vcf} | bgzip -c  > {job.output["vcf.bgz"]}',
#         )
#         job.command(f'tabix {job.output["vcf.bgz"]}')
#         get_batch().write_output(
#             job.output,
#             str(self.expected_outputs(multicohort)['wham_filtered_vcf']).replace('.vcf.bgz', ''),
#         )
#
#         expected_d = self.expected_outputs(multicohort)
#         return self.make_outputs(multicohort, data=expected_d, jobs=job)
#
#
# @stage(required_stages=[FilterWham], analysis_type='sv', analysis_keys=['annotated_vcf'])
# class AnnotateVcf(MultiCohortStage):
#     """
#     Add annotations, such as the inferred function and allele frequencies of variants,
#     to final VCF.
#
#     Annotations methods include:
#     * Functional annotation - annotate SVs with inferred functional consequence on
#       protein-coding regions, regulatory regions such as UTR and promoters, and other
#       non-coding elements.
#     * Allele frequency annotation - annotate SVs with their allele frequencies across
#       all samples, and samples of specific sex, as well as specific subpopulations.
#     * Allele Frequency annotation with external callset - annotate SVs with the allele
#       frequencies of their overlapping SVs in another callset, e.g. gnomad SV callset.
#
#     Note: the annotation stage is stupid, and re-orders the VCF by ID instead of position
#     This means that we run SVConcordance before this stage, to ensure that the IDs are correct
#     But we don't apply those IDs, in case multiple variants map to the same ID, which would
#     cause the variants to become unsorted when ordered by ID.
#     """
#
#     def expected_outputs(self, multicohort: MultiCohort) -> dict:
#         return {
#             'annotated_vcf': self.prefix / 'filtered_annotated.vcf.bgz',
#             'annotated_vcf_index': self.prefix / 'filtered_annotated.vcf.bgz.tbi',
#         }
#
#     def queue_jobs(self, multicohort: MultiCohort, inputs: StageInput) -> StageOutput | None:
#         """
#         Configure and queue jobs for SV annotation. Passing the VCF Index has become implicit
#         """
#         input_vcf = inputs.as_dict(multicohort, FilterWham)['wham_filtered_vcf']
#         expected_out = self.expected_outputs(multicohort)
#         billing_labels = {'stage': self.name.lower(), AR_GUID_NAME: try_get_ar_guid()}
#         job_or_none = queue_annotate_sv_jobs(multicohort, self.prefix, input_vcf, expected_out, billing_labels)
#         return self.make_outputs(multicohort, data=expected_out, jobs=job_or_none)
#
#
# @stage(required_stages=AnnotateVcf)
# class AnnotateVcfWithStrvctvre(MultiCohortStage):
#     def expected_outputs(self, multicohort: MultiCohort) -> dict[str, Path]:
#         return {
#             'strvctvre_vcf': self.prefix / 'strvctvre_annotated.vcf.gz',
#             'strvctvre_vcf_index': self.prefix / 'strvctvre_annotated.vcf.gz.tbi',
#         }
#
#     def queue_jobs(self, multicohort: MultiCohort, inputs: StageInput) -> StageOutput | None:
#
#         input_dict = inputs.as_dict(multicohort, AnnotateVcf)
#         expected_d = self.expected_outputs(multicohort)
#
#         # read vcf and index into the batch
#         input_vcf = get_batch().read_input_group(
#             vcf=str(input_dict['annotated_vcf']),
#             vcf_index=str(input_dict['annotated_vcf_index']),
#         )['vcf']
#
#         strvctvre_job = queue_annotate_strvctvre_job(input_vcf, str(expected_d['strvctvre_vcf']), self.get_job_attrs())
#
#         return self.make_outputs(multicohort, data=expected_d, jobs=strvctvre_job)
#
#
# @stage(required_stages=AnnotateVcfWithStrvctvre, analysis_type='sv')
# class SpiceUpSVIDs(MultiCohortStage):
#     """
#     Overwrites the GATK-SV assigned IDs with a meaningful ID
#     This new ID is either taken from an equivalent variant ID in the previous callset (found through SVConcordance)
#     or a new one is generated based on the call attributes itself
#     A boolean flag is used to switch behaviour, based on the truthiness of the return value of query_for_spicy_vcf
#     If a VCF is returned, True, which also means a VCF should have been used in UpdateStructuralVariantIDs
#     If the return is None, False, and the IDs will be generated from the call attributes alone
#     """
#
#     def expected_outputs(self, multicohort: MultiCohort) -> Path:
#         return self.prefix / 'fresh_ids.vcf.bgz'
#
#     def queue_jobs(self, multicohort: MultiCohort, inputs: StageInput) -> StageOutput:
#         input_vcf = get_batch().read_input(str(inputs.as_path(multicohort, AnnotateVcfWithStrvctvre, 'strvctvre_vcf')))
#         expected_output = self.expected_outputs(multicohort)
#
#         # update the IDs using a PythonJob
#         pyjob = get_batch().new_python_job('rename_sv_ids')
#         pyjob.storage('10Gi')
#         skip_prior_names = bool(query_for_spicy_vcf(multicohort.analysis_dataset.name))
#         pyjob.call(rename_sv_ids, input_vcf, pyjob.output, skip_prior_names)
#
#         # then compress & run tabix on that plain text result
#         bcftools_job = get_batch().new_job('bgzip and tabix')
#         bcftools_job.image(image_path('bcftools'))
#         bcftools_job.declare_resource_group(output={'vcf.bgz': '{root}.vcf.bgz', 'vcf.bgz.tbi': '{root}.vcf.bgz.tbi'})
#         bcftools_job.command(f'bcftools view {pyjob.output} | bgzip -c > {bcftools_job.output["vcf.bgz"]}')
#         bcftools_job.command(f'tabix {bcftools_job.output["vcf.bgz"]}')  # type: ignore
#
#         # get the output root to write to
#         get_batch().write_output(bcftools_job.output, str(expected_output).removesuffix('.vcf.bgz'))
#
#         return self.make_outputs(multicohort, data=expected_output, jobs=[pyjob, bcftools_job])
#
#
# @stage(required_stages=SpiceUpSVIDs)
# class AnnotateCohortSv(MultiCohortStage):
#     """
#     First step to transform annotated SV callset data into a seqr ready format
#     """
#
#     def expected_outputs(self, multicohort: MultiCohort) -> dict:
#         """
#         Expected to write a matrix table.
#         """
#         return {'tmp_prefix': str(self.tmp_prefix), 'mt': self.prefix / 'cohort_sv.mt'}
#
#     def queue_jobs(self, multicohort: MultiCohort, inputs: StageInput) -> StageOutput | None:
#         """
#         queue job(s) to rearrange the annotations prior to Seqr transformation
#         """
#         outputs = self.expected_outputs(multicohort)
#
#         vcf_path = inputs.as_path(target=multicohort, stage=SpiceUpSVIDs)
#         checkpoint_prefix = to_path(outputs['tmp_prefix']) / 'checkpoints'
#
#         job = annotate_cohort_jobs_sv(
#             vcf_path=vcf_path,
#             out_mt_path=outputs['mt'],
#             checkpoint_prefix=checkpoint_prefix,
#             job_attrs=self.get_job_attrs(multicohort),
#         )
#
#         return self.make_outputs(multicohort, data=outputs, jobs=job)
#
#
# @stage(required_stages=[CombineExclusionLists, AnnotateCohortSv], analysis_type='sv', analysis_keys=['mt'])
# class AnnotateDatasetSv(DatasetStage):
#     """
#     Subset the MT to be this Dataset only
#     Then work up all the genotype values
#     """
#
#     def expected_outputs(self, dataset: Dataset) -> dict[str, Path]:
#         """
#         Expected to generate a matrix table
#         """
#         return {'mt': (dataset.prefix() / 'mt' / f'SV-{get_workflow().output_version}-{dataset.name}.mt')}
#
#     def queue_jobs(self, dataset: Dataset, inputs: StageInput) -> StageOutput | None:
#         """
#         Subsets the whole MT to this cohort only
#         Then brings a range of genotype data into row annotations
#
#         Args:
#             dataset (Dataset): SGIDs specific to this dataset/project
#             inputs ():
#         """
#         # only create dataset MTs for datasets specified in the config
#         eligible_datasets = config_retrieve(['workflow', 'write_mt_for_datasets'], default=[])
#         if dataset.name not in eligible_datasets:
#             get_logger().info(f'Skipping AnnotateDatasetSv mt subsetting for {dataset}')
#             return None
#
#         mt_path = inputs.as_path(target=get_multicohort(), stage=AnnotateCohortSv, key='mt')
#         exclusion_file = inputs.as_path(get_multicohort(), stage=CombineExclusionLists, key='exclusion_list')
#
#         outputs = self.expected_outputs(dataset)
#
#         jobs = annotate_dataset_jobs_sv(
#             mt_path=mt_path,
#             sgids=dataset.get_sequencing_group_ids(),
#             out_mt_path=outputs['mt'],
#             tmp_prefix=self.tmp_prefix / dataset.name / 'checkpoints',
#             job_attrs=self.get_job_attrs(dataset),
#             exclusion_file=str(exclusion_file),
#         )
#
#         return self.make_outputs(dataset, data=outputs, jobs=jobs)
#
#
# @stage(
#     required_stages=[AnnotateDatasetSv],
#     analysis_type='es-index',  # specific type of es index
#     analysis_keys=['index_name'],
#     update_analysis_meta=lambda x: {'seqr-dataset-type': 'SV'},
# )
# class MtToEsSv(DatasetStage):
#     """
#     Create a Seqr index
#     https://github.com/populationgenomics/metamist/issues/539
#     """
#
#     def expected_outputs(self, dataset: Dataset) -> dict[str, str | Path]:
#         """
#         Expected to generate a Seqr index, which is not a file
#         """
#         sequencing_type = config_retrieve(['workflow', 'sequencing_type'])
#         index_name = f'{dataset.name}-{sequencing_type}-SV-{get_workflow().run_timestamp}'.lower()
#         return {
#             'index_name': index_name,
#             'done_flag': dataset.prefix() / 'es' / f'{index_name}.done',
#         }
#
#     def queue_jobs(self, dataset: Dataset, inputs: StageInput) -> StageOutput | None:
#         """
#         Uses the non-DataProc MT-to-ES conversion script
#         """
#         # only create the elasticsearch index for the datasets specified in the config
#         eligible_datasets = config_retrieve(['workflow', 'create_es_index_for_datasets'], default=[])
#         if dataset.name not in eligible_datasets:
#             get_logger().info(f'Skipping SV ES index creation for {dataset}')
#             return None
#
#         # try to generate a password here - we'll find out inside the script anyway, but
#         # by that point we'd already have localised the MT, wasting time and money
#         try:
#             _es_password_string = es_password()
#         except PermissionDenied:
#             get_logger().warning(f'No permission to access ES password, skipping for {dataset}')
#             return self.make_outputs(dataset)
#         except KeyError:
#             get_logger().warning(f'ES section not in config, skipping for {dataset}')
#             return self.make_outputs(dataset)
#
#         outputs = self.expected_outputs(dataset)
#
#         # get the absolute path to the MT
#         mt_path = str(inputs.as_path(target=dataset, stage=AnnotateDatasetSv, key='mt'))
#         # and just the name, used after localisation
#         mt_name = mt_path.split('/')[-1]
#
#         # get the expected outputs as Strings
#         index_name = str(outputs['index_name'])
#         flag_name = str(outputs['done_flag'])
#
#         job = get_batch().new_job(f'Generate {index_name} from {mt_path}')
#         # set all job attributes in one bash
#         job.cpu(4).memory('lowmem').storage('10Gi').image(config_retrieve(['workflow', 'driver_image']))
#
#         # localise the MT
#         job.command(f'gcloud --no-user-output-enabled storage cp -r {mt_path} $BATCH_TMPDIR')
#
#         # run the export from the localised MT - this job writes no new data, just transforms and exports over network
#         job.command(f'mt_to_es --mt_path "${{BATCH_TMPDIR}}/{mt_name}" --index {index_name} --flag {flag_name}')
#
#         return self.make_outputs(dataset, data=outputs['index_name'], jobs=job)
#
#
# @stage(required_stages=SpiceUpSVIDs, analysis_type='single_dataset_sv_annotated')
# class SplitAnnotatedSvVcfByDataset(DatasetStage):
#     """
#     takes the whole MultiCohort annotated VCF
#     splits it up into separate VCFs for each dataset
#     """
#
#     def expected_outputs(self, dataset: Dataset) -> Path:
#         return dataset.prefix() / 'annotated_sv.vcf.bgz'
#
#     def queue_jobs(self, dataset: Dataset, inputs: StageInput) -> StageOutput:
#         output = self.expected_outputs(dataset)
#         input_vcf = get_batch().read_input(inputs.as_path(get_multicohort(), SpiceUpSVIDs))
#
#         # write a temp file for this dataset containing all relevant SGs
#         sgids_list_path = dataset.tmp_prefix() / get_workflow().output_version / 'sgid-list.txt'
#         if not config_retrieve(['workflow', 'dry_run'], False):
#             with sgids_list_path.open('w') as f:
#                 for sgid in dataset.get_sequencing_group_ids():
#                     f.write(f'{sgid}\n')
#         job = get_batch().new_job(
#             name=f'SplitAnnotatedSvVcfByDataset: {dataset}',
#             attributes=self.get_job_attrs() | {'tool': 'bcftools'},
#         )
#         job.image(image_path('bcftools_120'))
#         job.cpu(1).memory('highmem').storage('10Gi')
#         job.declare_resource_group(output={'vcf.bgz': '{root}.vcf.bgz', 'vcf.bgz.tbi': '{root}.vcf.bgz.tbi'})
#
#         local_sgid_file = get_batch().read_input(sgids_list_path)
#         job.command(
#             f'bcftools view {input_vcf} '
#             f'--force-samples '
#             f'-S {local_sgid_file} '
#             f'-Oz -o {job.output["vcf.bgz"]} '
#             f'--write-index=tbi',
#         )
#
#         get_batch().write_output(job.output, str(output).removesuffix('.vcf.bgz'))
#
#         return self.make_outputs(dataset, data=output, jobs=job)


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
