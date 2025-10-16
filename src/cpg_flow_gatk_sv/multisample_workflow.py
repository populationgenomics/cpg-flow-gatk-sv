"""
All post-batching stages of the GATK-SV workflow
"""

import argparse

import loguru

from cpg_flow import stage, targets, workflow
from cpg_flow_gatk_sv import utils
from cpg_flow_gatk_sv.jobs.AnnotateCohort import create_annotate_cohort_job
from cpg_flow_gatk_sv.jobs.AnnotateDataset import create_annotate_dataset_jobs
from cpg_flow_gatk_sv.jobs.AnnotateVcf import create_svannotate_jobs
from cpg_flow_gatk_sv.jobs.AnnotateWithStrvctvre import create_strvctvre_jobs
from cpg_flow_gatk_sv.jobs.ClusterBatch import create_cluster_batch_jobs
from cpg_flow_gatk_sv.jobs.CombineExclusionLists import create_combine_exclusion_lists_job
from cpg_flow_gatk_sv.jobs.FilterBatch import create_filterbatch_jobs
from cpg_flow_gatk_sv.jobs.FilterGenotypes import create_filtergenotypes_jobs
from cpg_flow_gatk_sv.jobs.FilterWham import create_filter_wham_jobs
from cpg_flow_gatk_sv.jobs.FormatVcfForGatk import create_formatvcf_jobs
from cpg_flow_gatk_sv.jobs.GatherBatchEvidence import submit_gatherbatchevidence_jobs
from cpg_flow_gatk_sv.jobs.GenerateBatchMetrics import create_generate_batch_metrics_jobs
from cpg_flow_gatk_sv.jobs.GeneratePloidyTable import create_generate_ploidy_jobs
from cpg_flow_gatk_sv.jobs.GenotypeBatch import create_genotypebatch_jobs
from cpg_flow_gatk_sv.jobs.JoinRawCalls import create_joinrawcalls_jobs
from cpg_flow_gatk_sv.jobs.MakeCohortVcf import create_makecohortvcf_jobs
from cpg_flow_gatk_sv.jobs.MergeBatchSites import create_mergebatchsites_jobs
from cpg_flow_gatk_sv.jobs.MtToEs import create_mt_to_es_job
from cpg_flow_gatk_sv.jobs.SpiceUpSvIds import create_spicy_jobs
from cpg_flow_gatk_sv.jobs.SplitAnnotatedSvVcf import create_split_vcf_by_dataset_job
from cpg_flow_gatk_sv.jobs.SVConcordance import create_sv_concordance_jobs
from cpg_flow_gatk_sv.jobs.TrainGCNV import add_train_gcnv_jobs
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
class TrainGcnv(stage.CohortStage):
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
                for idx in range(config.config_retrieve(['workflow', 'model_tar_count']))
            ],
        }

    def queue_jobs(self, cohort: targets.Cohort, inputs: stage.StageInput) -> stage.StageOutput:
        outputs = self.expected_outputs(cohort)

        jobs = add_train_gcnv_jobs(cohort, output_dict=outputs)

        return self.make_outputs(cohort, data=outputs, jobs=jobs)


@stage.stage(required_stages=[MakeCohortCombinedPed, TrainGcnv])
class GatherBatchEvidence(stage.CohortStage):
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
        train_gcnv_outputs = inputs.as_dict(target=cohort, stage=TrainGcnv)

        outputs = self.expected_outputs(cohort)

        jobs = submit_gatherbatchevidence_jobs(
            pedigree_input=inputs.as_path(target=cohort, stage=MakeCohortCombinedPed),
            cohort=cohort,
            train_gcnv=train_gcnv_outputs,
            outputs=outputs,
        )

        return self.make_outputs(cohort, data=outputs, jobs=jobs)


@stage.stage(required_stages=[MakeCohortCombinedPed, GatherBatchEvidence])
class ClusterBatch(stage.CohortStage):
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
        batch_evidence_outputs = inputs.as_dict(cohort, GatherBatchEvidence)
        pedigree_input = inputs.as_str(target=cohort, stage=MakeCohortCombinedPed)

        outputs = self.expected_outputs(cohort)

        jobs = create_cluster_batch_jobs(
            pedigree_input=pedigree_input,
            cohort=cohort,
            batch_evidence_outputs=batch_evidence_outputs,
            outputs=outputs,
        )
        return self.make_outputs(cohort, data=outputs, jobs=jobs)


@stage.stage(required_stages=[MakeCohortCombinedPed, ClusterBatch, GatherBatchEvidence])
class GenerateBatchMetrics(stage.CohortStage):
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
        clusterbatch_outputs = inputs.as_dict(cohort, ClusterBatch)
        gatherbatchevidence_outputs = inputs.as_dict(cohort, GatherBatchEvidence)
        pedigree_input = inputs.as_str(target=cohort, stage=MakeCohortCombinedPed)

        outputs = self.expected_outputs(cohort)

        jobs = create_generate_batch_metrics_jobs(
            pedigree_input=pedigree_input,
            cohort=cohort,
            gather_batch_evidence_outputs=gatherbatchevidence_outputs,
            cluster_batch_outputs=clusterbatch_outputs,
            outputs=outputs,
        )

        return self.make_outputs(cohort, data=outputs, jobs=jobs)


@stage.stage(required_stages=[MakeCohortCombinedPed, GenerateBatchMetrics, ClusterBatch])
class FilterBatch(stage.CohortStage):
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
        metrics_d = inputs.as_dict(cohort, GenerateBatchMetrics)
        clusterbatch_d = inputs.as_dict(cohort, ClusterBatch)
        pedigree_input = inputs.as_str(target=cohort, stage=MakeCohortCombinedPed)

        outputs = self.expected_outputs(cohort)

        jobs = create_filterbatch_jobs(
            pedigree_input=pedigree_input,
            cohort=cohort,
            batch_metrics_outputs=metrics_d,
            cluster_batch_outputs=clusterbatch_d,
            outputs=outputs,
        )

        return self.make_outputs(cohort, data=outputs, jobs=jobs)


@stage.stage(required_stages=FilterBatch)
class MergeBatchSites(stage.MultiCohortStage):
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
        filter_batch_outputs = inputs.as_dict_by_target(FilterBatch)

        pesr_vcfs = [filter_batch_outputs[cohort.id]['filtered_pesr_vcf'] for cohort in multicohort.get_cohorts()]
        depth_vcfs = [filter_batch_outputs[cohort.id]['filtered_depth_vcf'] for cohort in multicohort.get_cohorts()]

        outputs = self.expected_outputs(multicohort)

        jobs = create_mergebatchsites_jobs(
            multicohort=multicohort,
            pesr_vcfs=pesr_vcfs,
            depth_vcfs=depth_vcfs,
            outputs=outputs,
        )

        return self.make_outputs(multicohort, data=outputs, jobs=jobs)


@stage.stage(analysis_type='sv', required_stages=FilterBatch)
class CombineExclusionLists(stage.MultiCohortStage):
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

        filter_batch_outputs = inputs.as_dict_by_target(FilterBatch)
        all_filter_lists = [
            str(filter_batch_outputs[cohort.id]['outlier_samples_excluded_file'])
            for cohort in multicohort.get_cohorts()
        ]

        output = self.expected_outputs(multicohort)

        job = create_combine_exclusion_lists_job(all_filter_lists, output)

        return self.make_outputs(multicohort, data=output, jobs=job)


@stage.stage(
    required_stages=[FilterBatch, GatherBatchEvidence, MergeBatchSites],
    analysis_type='sv',
    analysis_keys=['genotyped_depth_vcf', 'genotyped_pesr_vcf'],
)
class GenotypeBatch(stage.CohortStage):
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
        filterbatch_outputs = inputs.as_dict(cohort, FilterBatch)
        batchevidence_outputs = inputs.as_dict(cohort, GatherBatchEvidence)
        mergebatch_outputs = inputs.as_dict(workflow.get_multicohort(), MergeBatchSites)

        outputs = self.expected_outputs(cohort)

        jobs = create_genotypebatch_jobs(
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
        GatherBatchEvidence,
        GenotypeBatch,
        FilterBatch,
    ]
)
class MakeCohortVcf(stage.MultiCohortStage):
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
        gatherbatchevidence_outputs = inputs.as_dict_by_target(GatherBatchEvidence)
        genotypebatch_outputs = inputs.as_dict_by_target(GenotypeBatch)
        filterbatch_outputs = inputs.as_dict_by_target(FilterBatch)
        pedigree_input = inputs.as_path(target=multicohort, stage=MakeMultiCohortCombinedPed)

        outputs = self.expected_outputs(multicohort)

        jobs = create_makecohortvcf_jobs(
            multicohort,
            pedigree_input,
            gatherbatchevidence_outputs,
            genotypebatch_outputs,
            filterbatch_outputs,
            outputs,
        )

        return self.make_outputs(multicohort, data=outputs, jobs=jobs)


@stage.stage(required_stages=[MakeMultiCohortCombinedPed, MakeCohortVcf])
class FormatVcfForGatk(stage.MultiCohortStage):
    def expected_outputs(self, multicohort: targets.MultiCohort) -> dict:
        return {
            'gatk_formatted_vcf': self.prefix / 'gatk_formatted.vcf.gz',
            'gatk_formatted_vcf_index': self.prefix / 'gatk_formatted.vcf.gz.tbi',
        }

    def queue_jobs(self, multicohort: targets.MultiCohort, inputs: stage.StageInput) -> stage.StageOutput:
        pedigree_input = inputs.as_str(target=multicohort, stage=MakeMultiCohortCombinedPed)
        vcf_input = inputs.as_str(multicohort, MakeCohortVcf, 'vcf')
        outputs = self.expected_outputs(multicohort)

        jobs = create_formatvcf_jobs(
            multicohort,
            pedigree_input,
            vcf_input,
            outputs,
        )

        return self.make_outputs(multicohort, data=outputs, jobs=jobs)


@stage.stage(required_stages=[MakeMultiCohortCombinedPed, ClusterBatch])
class JoinRawCalls(stage.MultiCohortStage):
    """
    Joins all individually clustered caller results
    """

    def expected_outputs(self, multicohort: targets.MultiCohort) -> dict[str, Path]:
        return {
            'joined_raw_calls_vcf': self.prefix / 'raw_clustered_calls.vcf.gz',
            'joined_raw_calls_vcf_index': self.prefix / 'raw_clustered_calls.vcf.gz.tbi',
        }

    def queue_jobs(self, multicohort: targets.MultiCohort, inputs: stage.StageInput) -> stage.StageOutput:
        pedigree = inputs.as_str(target=multicohort, stage=MakeMultiCohortCombinedPed)
        clusterbatch_outputs = inputs.as_dict_by_target(ClusterBatch)

        outputs = self.expected_outputs(multicohort)

        jobs = create_joinrawcalls_jobs(
            multicohort=multicohort,
            pedigree=pedigree,
            clusterbatch_outputs=clusterbatch_outputs,
            outputs=outputs,
        )

        return self.make_outputs(multicohort, data=outputs, jobs=jobs)


@stage.stage(required_stages=[JoinRawCalls, FormatVcfForGatk])
class SVConcordance(stage.MultiCohortStage):
    """
    Takes the clean VCF and reformat for GATK intake
    """

    def expected_outputs(self, multicohort: targets.MultiCohort) -> dict:
        """
        create dictionary of names -> output paths
        """

        return {
            'concordance_vcf': self.prefix / 'sv_concordance.vcf.gz',
            'concordance_vcf_index': self.prefix / 'sv_concordance.vcf.gz.tbi',
        }

    def queue_jobs(self, multicohort: targets.MultiCohort, inputs: stage.StageInput) -> stage.StageOutput:
        """
        configure and queue jobs for SV concordance
        """
        outputs = self.expected_outputs(multicohort)

        formatvcf_output = inputs.as_str(multicohort, FormatVcfForGatk, 'gatk_formatted_vcf')
        joinrawcalls_output = inputs.as_str(multicohort, JoinRawCalls, 'joined_raw_calls_vcf')

        jobs = create_sv_concordance_jobs(
            multicohort=multicohort,
            formatvcf_output=formatvcf_output,
            joinrawcalls_output=joinrawcalls_output,
            outputs=outputs,
        )

        return self.make_outputs(multicohort, data=outputs, jobs=jobs)


@stage.stage(required_stages=[MakeMultiCohortCombinedPed, SVConcordance])
class GeneratePloidyTable(stage.MultiCohortStage):
    """
    Quick PythonJob to generate a ploidy table
    Calls a homebrewed version of this table generator:
    github.com/broadinstitute/gatk-sv/blob/main/src/sv-pipeline/scripts/ploidy_table_from_ped.py
    """

    def expected_outputs(self, multicohort: targets.MultiCohort) -> Path:
        """
        only one output, the ploidy table
        """

        return self.prefix / 'ploidy_table.txt'

    def queue_jobs(self, multicohort: targets.MultiCohort, inputs: stage.StageInput) -> stage.StageOutput:
        output = self.expected_outputs(multicohort)
        pedigree_input = inputs.as_path(target=multicohort, stage=MakeMultiCohortCombinedPed)

        job = create_generate_ploidy_jobs(
            pedigree=pedigree_input,
            output=output,
        )

        return self.make_outputs(multicohort, data=output, jobs=job)


@stage.stage(
    required_stages=[MakeMultiCohortCombinedPed, GeneratePloidyTable, SVConcordance],
)
class FilterGenotypes(stage.MultiCohortStage):
    """
    Steps required to post-filter called genotypes
    """

    def expected_outputs(self, multicohort: targets.MultiCohort) -> dict:
        return {
            'filtered_vcf': self.prefix / 'filtered.vcf.gz',
            'filtered_vcf_index': self.prefix / 'filtered.vcf.gz.tbi',
            'unfiltered_recalibrated_vcf': self.prefix / 'unfiltered_recalibrated.vcf.gz',
            'unfiltered_recalibrated_vcf_index': self.prefix / 'unfiltered_recalibrated.vcf.gz.tbi',
        }

    def queue_jobs(self, multicohort: targets.MultiCohort, inputs: stage.StageInput) -> stage.StageOutput:
        outputs = self.expected_outputs(multicohort)

        pedigree_input = inputs.as_str(target=multicohort, stage=MakeMultiCohortCombinedPed)
        ploidy_table = inputs.as_str(multicohort, GeneratePloidyTable)
        sv_conc_vcf = inputs.as_str(multicohort, SVConcordance, 'concordance_vcf')

        jobs = create_filtergenotypes_jobs(
            multicohort=multicohort,
            pedigree=pedigree_input,
            ploidy_table=ploidy_table,
            sv_conc_vcf=sv_conc_vcf,
            outputs=outputs,
        )

        return self.make_outputs(multicohort, data=outputs, jobs=jobs)


@stage.stage(required_stages=FilterGenotypes)
class UpdateStructuralVariantIDs(stage.MultiCohortStage):
    """
    Runs SVConcordance between the results of this callset and the results of a previous callset
    This causes the Variant IDs of a matching variant to be updated to the previous callset's ID
    Consistency of Variant ID is crucial to Seqr/AIP identifying the same variant across different callsets
    By default GATK-SV creates an auto-incrementing ID, rather than one based on variant attributes
    If a new call is added at a chromosome start for any sample, the ID for all variants on the chromosome
    will be shifted up, meaning that Seqr labels/any other process linked to the ID will be incorrect
    """

    def expected_outputs(self, multicohort: targets.MultiCohort) -> dict[str, Path]:
        return {
            'concordance_vcf': self.prefix / 'updated_ids.vcf.gz',
            'concordance_vcf_index': self.prefix / 'updated_ids.vcf.gz.tbi',
        }

    def queue_jobs(self, multicohort: targets.MultiCohort, inputs: stage.StageInput) -> stage.StageOutput:
        # allow for no prior name/IDs
        if not (spicy_vcf := utils.query_for_spicy_vcf(multicohort.analysis_dataset.name)):
            loguru.logger.info('No previous Spicy VCF found for {cohort.analysis_dataset.name}')
            return self.make_outputs(multicohort, skipped=True)

        outputs = self.expected_outputs(multicohort)

        filter_genotypes_output = inputs.as_str(multicohort, FilterGenotypes, key='filtered_vcf')

        jobs = create_sv_concordance_jobs(
            multicohort=multicohort,
            formatvcf_output=filter_genotypes_output,
            joinrawcalls_output=spicy_vcf,
            outputs=outputs,
        )

        return self.make_outputs(multicohort, data=outputs, jobs=jobs)


@stage.stage(
    required_stages=[FilterGenotypes, UpdateStructuralVariantIDs],
    analysis_type='sv',
)
class FilterWham(stage.MultiCohortStage):
    """
    Filters the VCF to remove deletions only called by Wham
    github.com/broadinstitute/gatk-sv/blob/main/wdl/ApplyManualVariantFilter.wdl
    """

    def expected_outputs(self, multicohort: targets.MultiCohort) -> Path:
        return self.prefix / 'filtered.vcf.bgz'

    def queue_jobs(self, multicohort: targets.MultiCohort, inputs: stage.StageInput) -> stage.StageOutput:
        """
        configure and queue jobs for SV annotation
        passing the VCF Index has become implicit, which may be a problem for us
        """
        # read the concordance-with-prev-batch VCF if appropriate, otherwise use the filtered VCF
        if utils.query_for_spicy_vcf(multicohort.analysis_dataset.name):
            loguru.logger.info(f'Variant IDs were updated for {multicohort.analysis_dataset.name}')
            input_vcf = inputs.as_str(multicohort, UpdateStructuralVariantIDs, 'concordance_vcf')
        else:
            loguru.logger.info(f'No Spicy VCF was found, default IDs for {multicohort.analysis_dataset.name}')
            input_vcf = inputs.as_str(multicohort, FilterGenotypes, 'filtered_vcf')

        output = self.expected_outputs(multicohort)

        job = create_filter_wham_jobs(input_vcf, output=str(output))

        return self.make_outputs(multicohort, data=output, jobs=job)


@stage.stage(
    required_stages=[FilterWham, MakeMultiCohortCombinedPed],
    analysis_type='sv',
    analysis_keys=['annotated_vcf'],
)
class AnnotateVcf(stage.MultiCohortStage):
    """
    Add annotations, such as the inferred function and allele frequencies of variants,
    to final VCF.

    Annotations methods include:
    * Functional annotation - annotate SVs with inferred functional consequence on
      protein-coding regions, regulatory regions such as UTR and promoters, and other
      non-coding elements.
    * Allele frequency annotation - annotate SVs with their allele frequencies across
      all samples, and samples of specific sex, as well as specific subpopulations.
    * Allele Frequency annotation with external callset - annotate SVs with the allele
      frequencies of their overlapping SVs in another callset, e.g. gnomad SV callset.

    Note: the annotation stage is stupid, and re-orders the VCF by ID instead of position
    This means that we run SVConcordance before this stage, to ensure that the IDs are correct
    But we don't apply those IDs, in case multiple variants map to the same ID, which would
    cause the variants to become unsorted when ordered by ID.
    """

    def expected_outputs(self, multicohort: targets.MultiCohort) -> dict[str, Path]:
        return {
            'annotated_vcf': self.prefix / 'filtered_annotated.vcf.bgz',
            'annotated_vcf_index': self.prefix / 'filtered_annotated.vcf.bgz.tbi',
        }

    def queue_jobs(self, multicohort: targets.MultiCohort, inputs: stage.StageInput) -> stage.StageOutput:
        """
        Configure and queue jobs for SV annotation. Passing the VCF Index has become implicit
        """
        outputs = self.expected_outputs(multicohort)

        input_vcf = inputs.as_str(multicohort, FilterWham)
        pedigree = inputs.as_str(target=multicohort, stage=MakeMultiCohortCombinedPed)

        jobs = create_svannotate_jobs(
            multicohort=multicohort,
            input_vcf=input_vcf,
            pedigree=pedigree,
            outputs=outputs,
        )

        return self.make_outputs(multicohort, data=outputs, jobs=jobs)


@stage.stage(required_stages=AnnotateVcf)
class AnnotateVcfWithStrvctvre(stage.MultiCohortStage):
    def expected_outputs(self, multicohort: targets.MultiCohort) -> dict[str, Path]:
        return {
            'strvctvre_vcf': self.prefix / 'strvctvre_annotated.vcf.gz',
            'strvctvre_vcf_index': self.prefix / 'strvctvre_annotated.vcf.gz.tbi',
        }

    def queue_jobs(self, multicohort: targets.MultiCohort, inputs: stage.StageInput) -> stage.StageOutput:
        input_vcf = inputs.as_str(multicohort, AnnotateVcf, 'annotated_vcf')
        outputs = self.expected_outputs(multicohort)

        job = create_strvctvre_jobs(
            input_vcf=input_vcf,
            output=str(outputs['strvctvre_vcf']),
            name=self.name,
        )
        return self.make_outputs(multicohort, data=outputs, jobs=job)


@stage.stage(required_stages=AnnotateVcfWithStrvctvre, analysis_type='sv')
class SpiceUpSvIds(stage.MultiCohortStage):
    """
    Overwrites the GATK-SV assigned IDs with a meaningful ID
    This new ID is either taken from an equivalent variant ID in the previous callset (found through SVConcordance)
    or a new one is generated based on the call attributes itself
    A boolean flag is used to switch behaviour, based on the truthiness of the return value of query_for_spicy_vcf
    If a VCF is returned, True, which also means a VCF should have been used in UpdateStructuralVariantIDs
    If the return is None, False, and the IDs will be generated from the call attributes alone
    """

    def expected_outputs(self, multicohort: targets.MultiCohort) -> Path:
        return self.prefix / 'fresh_ids.vcf.bgz'

    def queue_jobs(self, multicohort: targets.MultiCohort, inputs: stage.StageInput) -> stage.StageOutput:
        output = self.expected_outputs(multicohort)
        input_vcf = inputs.as_str(multicohort, AnnotateVcfWithStrvctvre, 'strvctvre_vcf')

        jobs = create_spicy_jobs(
            input_vcf=input_vcf,
            skip_prior_names=bool(utils.query_for_spicy_vcf(multicohort.analysis_dataset.name)),
            output=str(output),
        )

        return self.make_outputs(multicohort, data=output, jobs=jobs)


@stage.stage(required_stages=SpiceUpSvIds)
class AnnotateCohort(stage.MultiCohortStage):
    """
    First step to transform annotated SV callset data into a seqr ready format
    """

    def expected_outputs(self, multicohort: targets.MultiCohort) -> Path:
        """
        Expected to write a matrix table.
        """
        return self.prefix / 'cohort_sv.mt'

    def queue_jobs(self, multicohort: targets.MultiCohort, inputs: stage.StageInput) -> stage.StageOutput:
        """
        queue job(s) to rearrange the annotations prior to Seqr transformation
        """
        output = self.expected_outputs(multicohort)

        vcf_path = inputs.as_str(target=multicohort, stage=SpiceUpSvIds)

        job = create_annotate_cohort_job(
            vcf=vcf_path,
            out_mt=output,
            checkpoint=self.tmp_prefix / 'checkpoint.mt',
            job_attrs=self.get_job_attrs(multicohort),
        )

        return self.make_outputs(multicohort, data=output, jobs=job)


@stage.stage(
    required_stages=[CombineExclusionLists, AnnotateCohort],
    analysis_type='sv',
)
class AnnotateDataset(stage.DatasetStage):
    """
    Subset the MT to be this Dataset only
    Then work up all the genotype values
    """

    def expected_outputs(self, dataset: targets.Dataset) -> Path:
        return dataset.prefix() / 'mt' / f'SV-{workflow.get_workflow().output_version}-{dataset.name}.mt'

    def queue_jobs(self, dataset: targets.Dataset, inputs: stage.StageInput) -> stage.StageOutput | None:
        """
        Subsets the whole MT to this cohort only
        Then brings a range of genotype data into row annotations

        Args:
            dataset (Dataset): SGIDs specific to this dataset/project
            inputs ():
        """
        # only create dataset MTs for datasets specified in the config
        if dataset.name not in config.config_retrieve(['workflow', 'write_mt_for_datasets'], default=[]):
            loguru.logger.info(f'Skipping AnnotateDatasetStage mt subsetting for {dataset.name}')
            return None

        multicohort = workflow.get_multicohort()
        cohort_mt = inputs.as_str(target=multicohort, stage=AnnotateCohort)
        exclusion_file = inputs.as_str(multicohort, stage=CombineExclusionLists)

        output = self.expected_outputs(dataset)

        dataset_sgid_file = utils.write_dataset_sg_ids(dataset)

        job = create_annotate_dataset_jobs(
            mt=cohort_mt,
            sgid_file=dataset_sgid_file,
            mt_out=output,
            dataset_mt=self.tmp_prefix / f'{dataset.name}_subset.mt',
            job_attrs=self.get_job_attrs(dataset),
            exclusion_file=exclusion_file,
        )

        return self.make_outputs(dataset, data=output, jobs=job)


@stage.stage(
    required_stages=[AnnotateDataset],
    analysis_type='es-index',
    analysis_keys=['index_name'],
    update_analysis_meta=lambda x: {'seqr-dataset-type': 'SV'},  # noqa: ARG005
)
class MtToEs(stage.DatasetStage):
    """
    Create a Seqr index
    https://github.com/populationgenomics/metamist/issues/539
    """

    def expected_outputs(self, dataset: targets.Dataset) -> dict[str, str | Path]:
        """
        Expected to generate a Seqr index, which is not a file
        """
        sequencing_type = config.config_retrieve(['workflow', 'sequencing_type'])
        index_name = f'{dataset.name}-{sequencing_type}-SV-{workflow.get_workflow().run_timestamp}'.lower()
        return {
            'index_name': index_name,
            'done_flag': dataset.prefix() / 'es' / f'{index_name}.done',
        }

    def queue_jobs(self, dataset: targets.Dataset, inputs: stage.StageInput) -> stage.StageOutput | None:
        """
        Uses the non-DataProc MT-to-ES conversion script
        """
        # only create the elasticsearch index for the datasets specified in the config
        if dataset.name not in config.config_retrieve(['workflow', 'create_es_index_for_datasets'], default=[]):
            loguru.logger.info(f'Skipping SV ES index creation for {dataset}')
            return None

        outputs = self.expected_outputs(dataset)

        # get the absolute path to the MT
        mt_path = inputs.as_str(target=dataset, stage=AnnotateDataset)
        job_or_none = create_mt_to_es_job(
            dataset=dataset,
            mt_path=mt_path,
            outputs=outputs,
            job_attrs=self.get_job_attrs(target=dataset) | {'tool': 'hail'},
        )

        return self.make_outputs(dataset, data=outputs, jobs=job_or_none)


@stage.stage(
    required_stages=SpiceUpSvIds,
    analysis_type='single_dataset_sv_annotated',
)
class SplitAnnotatedSvVcfByDataset(stage.DatasetStage):
    """
    takes the whole MultiCohort annotated VCF
    splits it up into separate VCFs for each dataset
    """

    def expected_outputs(self, dataset: targets.Dataset) -> Path:
        return (
            dataset.prefix() / 'gatk_sv' / workflow.get_workflow().output_version / self.name / 'annotated_sv.vcf.bgz'
        )

    def queue_jobs(self, dataset: targets.Dataset, inputs: stage.StageInput) -> stage.StageOutput:
        output = self.expected_outputs(dataset)

        input_vcf = inputs.as_str(workflow.get_multicohort(), SpiceUpSvIds)

        dataset_sgid_file = utils.write_dataset_sg_ids(dataset)

        job = create_split_vcf_by_dataset_job(
            dataset=dataset,
            input_vcf=input_vcf,
            dataset_sgid_file=dataset_sgid_file,
            output=str(output),
            job_attrs=self.get_job_attrs(dataset),
        )

        return self.make_outputs(dataset, data=output, jobs=job)


def cli_main():
    """
    CLI entrypoint - starts up the workflow
    """
    parser = argparse.ArgumentParser()
    parser.add_argument('--dry_run', action='store_true', help='Dry run')
    args = parser.parse_args()

    workflow.run_workflow(
        name='gatk_sv_multisample',
        stages=[
            TrainGcnv,
            MakeCohortCombinedPed,
            MakeMultiCohortCombinedPed,
            GatherBatchEvidence,
            ClusterBatch,
            GenerateBatchMetrics,
            FilterBatch,
            MergeBatchSites,
            CombineExclusionLists,
            GenotypeBatch,
            MakeCohortVcf,
            JoinRawCalls,
            GeneratePloidyTable,
            FilterGenotypes,
            AnnotateVcf,
            AnnotateVcfWithStrvctvre,
            AnnotateCohort,
            AnnotateDataset,
            MtToEs,
            SplitAnnotatedSvVcfByDataset,
        ],
        dry_run=args.dry_run,
    )


if __name__ == '__main__':
    cli_main()
