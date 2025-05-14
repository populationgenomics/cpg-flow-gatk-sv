"""
single-sample components of the GATK SV workflow
"""

from typing import TYPE_CHECKING
import json
import logging
from functools import cache

from cpg_utils.config import config_retrieve
from cpg_flow_gatk_sv.utils import (
    SV_CALLERS,
)
from cpg_flow.stage import CohortStage, MultiCohortStage, SequencingGroupStage, stage

from cpg_flow_gatk_sv.jobs.GatherSampleEvidence import create_gather_sample_evidence_jobs

if TYPE_CHECKING:
    from cpg_utils import Path
    from cpg_flow.targets import SequencingGroup
    from cpg_flow.stage import StageInput, StageOutput


@cache
def get_sv_callers():
    if only_jobs := config_retrieve(['workflow', 'GatherSampleEvidence', 'only_jobs'], None):
        callers = [caller for caller in SV_CALLERS if caller in only_jobs]
        if not callers:
            logging.warning('No SV callers enabled')
        return callers
    return SV_CALLERS


@stage(analysis_keys=[f'{caller}_vcf' for caller in get_sv_callers()] if get_sv_callers() else None, analysis_type='sv')
class GatherSampleEvidence(SequencingGroupStage):
    """
    https://github.com/broadinstitute/gatk-sv#gathersampleevidence
    https://github.com/broadinstitute/gatk-sv/blob/master/wdl/GatherSampleEvidence.wdl
    """

    def expected_outputs(self, sequencing_group: 'SequencingGroup') -> 'dict[str, Path]':
        """
        Expected to produce coverage counts, a VCF for each variant caller,
        and a txt for each type of SV evidence (SR, PE, SD).

        it's critical to separate the ending with a dot, e.g.: `*.sr.txt.gz`,
        These files are passed to `gatk PrintSVEvidence`, that determines file
        format based on the file name.
        It would strongly expect the files to end exactly with either
        `.sr.txt.gz`, `.pe.txt.gz`, or `.sd.txt.gz`, otherwise it would fail with
        "A USER ERROR has occurred: Cannot read file:///cromwell_root/... because
        no suitable codecs found".
        """

        prefix = sequencing_group.make_sv_evidence_path

        outputs: dict[str, 'Path'] = {
            'coverage_counts': prefix / f'{sequencing_group.id}.coverage_counts.tsv.gz',
            # split reads
            'pesr_split': prefix / f'{sequencing_group.id}.sr.txt.gz',
            'pesr_split_index': prefix / f'{sequencing_group.id}.sr.txt.gz.tbi',
            # site depth
            'pesr_sd': prefix / f'{sequencing_group.id}.sd.txt.gz',
            'pesr_sd_index': prefix / f'{sequencing_group.id}.sd.txt.gz.tbi',
            # discordant paired reads
            'pesr_disc': prefix / f'{sequencing_group.id}.pe.txt.gz',
            'pesr_disc_index': prefix / f'{sequencing_group.id}.pe.txt.gz.tbi',
        }

        # Caller's VCFs
        for caller in get_sv_callers():
            outputs[f'{caller}_vcf'] = prefix / f'{sequencing_group.id}.{caller}.vcf.gz'
            outputs[f'{caller}_index'] = prefix / f'{sequencing_group.id}.{caller}.vcf.gz.tbi'

        # TODO This selection process may need to adapt to a new condition...
        # TODO If Scramble is being run, but Manta is not, manta_vcf and index becomes a required input
        if only_jobs := config_retrieve(['workflow', self.name, 'only_jobs'], None):
            # remove the expected outputs for the jobs that are not in only_jobs
            new_expected = {}
            for job in only_jobs:
                for key, path in outputs.items():
                    if job in key:
                        new_expected[key] = path
            outputs = new_expected

        return outputs

    def queue_jobs(self, sequencing_group: 'SequencingGroup', inputs: 'StageInput') -> 'StageOutput':
        """
        Add jobs to batch
        Adds billing-related labels to the Cromwell job(s)
        """

        outputs = self.expected_outputs(sequencing_group)

        jobs = create_gather_sample_evidence_jobs(
            sg=sequencing_group,
            expected_outputs=outputs,
        )
        return self.make_outputs(sequencing_group, data=outputs, jobs=jobs)


# @stage(required_stages=GatherSampleEvidence, analysis_type='qc', analysis_keys=['qc_table'])
# class EvidenceQC(CohortStage):
#     """
#     https://github.com/broadinstitute/gatk-sv#evidenceqc
#     """
#
#     def expected_outputs(self, cohort: Cohort) -> dict[str, Path]:
#         """
#         Expected to return a bunch of batch-level summary files.
#         """
#         fname_by_key = {
#             'ploidy_matrix': 'ploidy_matrix.bed.gz',
#             'ploidy_plots': 'ploidy_plots.tar.gz',
#             'WGD_dist': 'WGD_score_distributions.pdf',
#             'WGD_matrix': 'WGD_scoring_matrix_output.bed.gz',
#             'WGD_scores': 'WGD_scores.txt.gz',
#             'bincov_matrix': f'{self.name}.RD.txt.gz',
#             'bincov_matrix_index': f'{self.name}.RD.txt.gz.tbi',
#             'bincov_median': 'medianCov.transposed.bed',
#             'qc_table': 'evidence_qc_table.tsv',
#         }
#         for caller in SV_CALLERS:
#             for k in ['low', 'high']:
#                 fname_by_key[f'{caller}_qc_{k}'] = f'{caller}_QC.outlier.{k}'
#
#         return {key: self.get_stage_cohort_prefix(cohort) / fname for key, fname in fname_by_key.items()}
#
#     def queue_jobs(self, cohort: Cohort, inputs: StageInput) -> StageOutput:
#         d = inputs.as_dict_by_target(GatherSampleEvidence)
#         sgids = cohort.get_sequencing_group_ids()
#
#         input_dict: dict[str, Any] = {
#             'batch': cohort.id,
#             'samples': sgids,
#             'run_vcf_qc': True,
#             'counts': [str(d[sid]['coverage_counts']) for sid in sgids],
#         }
#         for caller in SV_CALLERS:
#             input_dict[f'{caller}_vcfs'] = [str(d[sid][f'{caller}_vcf']) for sid in sgids]
#
#         input_dict |= get_images(
#             ['sv_base_mini_docker', 'sv_base_docker', 'sv_pipeline_docker', 'sv_pipeline_qc_docker'],
#         )
#
#         input_dict |= get_references(['genome_file', 'wgd_scoring_mask'])
#
#         expected_d = self.expected_outputs(cohort)
#
#         billing_labels = {'stage': self.name.lower(), AR_GUID_NAME: try_get_ar_guid()}
#
#         # runs for approx 5 hours, depending on sample count
#         jobs = add_gatk_sv_jobs(
#             dataset=cohort.analysis_dataset,
#             wfl_name=self.name,
#             input_dict=input_dict,
#             expected_out_dict=expected_d,
#             labels=billing_labels,
#             job_size=CromwellJobSizes.MEDIUM,
#         )
#         return self.make_outputs(cohort, data=expected_d, jobs=jobs)
#
#
# @stage(required_stages=EvidenceQC, analysis_type='sv', analysis_keys=['batch_json'])
# class CreateSampleBatches(MultiCohortStage):
#     """
#     uses the values generated in EvidenceQC
#     splits the sequencing groups into batches based on median coverage,
#     PCR +/- status, and Sex
#
#     The output of this Stage will contain the distinct SG batches to use for the
#     following series of Stages. For now, standard practice is to create a separate
#     minimal configuration file for each sub-batch, containing the list of SG IDs
#     as the `only_sgs` key. The gatk_sv_multisample_1 and gatk_sv_sandwich WFs are
#     then run separately for each sub-batch, with the active SGs controlled via the
#     config contents.
#
#     The output of this stage is used to generate custom cohorts
#     """
#
#     def expected_outputs(self, multicohort: MultiCohort) -> dict:
#         return {'batch_json': self.prefix / 'sgid_batches.json'}
#
#     def queue_jobs(self, multicohort: MultiCohort, inputs: StageInput) -> StageOutput | None:
#         """
#         this stage has been clarified - this will only run on Genomes
#         Exomes were never a supported case, they have a separate pipeline
#         """
#
#         expected = self.expected_outputs(multicohort)
#
#         if config_retrieve(['workflow', 'sequencing_type']) != 'genome':
#             raise RuntimeError('This workflow is not intended for Exome data')
#
#         # get the batch size parameters
#         min_batch_size = config_retrieve(['workflow', 'min_batch_size'], 100)
#         max_batch_size = config_retrieve(['workflow', 'max_batch_size'], 300)
#
#         # Get the sequencing groups
#         sequencing_groups = {
#             sequencing_group.id: sequencing_group.meta for sequencing_group in multicohort.get_sequencing_groups()
#         }
#         if len(sequencing_groups) < min_batch_size:
#             logging.error('Too few sequencing groups to form batches')
#             raise RuntimeError('too few samples to create batches')
#
#         # write them to a json file in tmp
#         sgs_json_path = to_path(self.tmp_prefix / 'sgs_meta.json')
#         with sgs_json_path.open('w') as f:
#             json.dump(sequencing_groups, f)
#
#         # Get the QC tables for each input cohort
#         qc_tables = [inputs.as_dict(cohort, EvidenceQC)['qc_table'] for cohort in multicohort.get_cohorts()]
#
#         py_job = get_batch().new_python_job('create_sample_batches')
#         py_job.image(config_retrieve(['workflow', 'driver_image']))
#         py_job.call(
#             sample_batching.partition_batches,
#             qc_tables,
#             sgs_json_path,
#             str(expected['batch_json']),
#             min_batch_size,
#             max_batch_size,
#         )
#
#         return self.make_outputs(multicohort, data=expected, jobs=py_job)
