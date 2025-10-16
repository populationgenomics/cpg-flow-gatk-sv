"""
single-sample components of the GATK SV workflow
"""

import argparse

from cpg_flow import stage, targets, workflow
from cpg_flow_gatk_sv import utils
from cpg_flow_gatk_sv.jobs.CreateSampleBatches import create_sample_batches
from cpg_flow_gatk_sv.jobs.EvidenceQC import create_evidence_qc_jobs
from cpg_flow_gatk_sv.jobs.GatherSampleEvidence import create_gather_sample_evidence_jobs
from cpg_utils import Path, config


@stage.stage(
    analysis_keys=[f'{caller}_vcf' for caller in utils.get_sv_callers()] if utils.get_sv_callers() else None,
    analysis_type='sv',
)
class GatherSampleEvidence(stage.SequencingGroupStage):
    """
    https://github.com/broadinstitute/gatk-sv#gathersampleevidence
    https://github.com/broadinstitute/gatk-sv/blob/master/wdl/GatherSampleEvidence.wdl
    """

    def expected_outputs(self, sequencing_group: targets.SequencingGroup) -> dict[str, Path]:
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

        outputs: dict[str, Path] = {
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
        for caller in utils.get_sv_callers(add_manta=True):
            outputs[f'{caller}_vcf'] = prefix / f'{sequencing_group.id}.{caller}.vcf.gz'
            outputs[f'{caller}_index'] = prefix / f'{sequencing_group.id}.{caller}.vcf.gz.tbi'

        if only_jobs := config.config_retrieve(['workflow', self.name, 'only_jobs'], None):
            # remove the expected outputs for the jobs that are not in only_jobs
            new_expected = {}
            for job in only_jobs:
                for key, path in outputs.items():
                    if job in key:
                        new_expected[key] = path

            if 'scramble' in only_jobs and 'manta' not in only_jobs:
                # if Scramble is being run, but Manta is not, manta_vcf and index becomes a required input
                # so we need to add them to the expected outputs
                new_expected['manta_vcf'] = outputs['manta_vcf']
                new_expected['manta_index'] = outputs['manta_index']

            outputs = new_expected

        return outputs

    def queue_jobs(self, sequencing_group: targets.SequencingGroup, inputs: stage.StageInput) -> stage.StageOutput:
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


@stage.stage(required_stages=GatherSampleEvidence, analysis_type='qc', analysis_keys=['qc_table'])
class EvidenceQC(stage.CohortStage):
    """
    https://github.com/broadinstitute/gatk-sv#evidenceqc
    """

    def expected_outputs(self, cohort: targets.Cohort) -> dict[str, Path]:
        """
        Expected to return a bunch of batch-level summary files.
        """
        fname_by_key = {
            'ploidy_matrix': 'ploidy_matrix.bed.gz',
            'ploidy_plots': 'ploidy_plots.tar.gz',
            'WGD_dist': 'WGD_score_distributions.pdf',
            'WGD_matrix': 'WGD_scoring_matrix_output.bed.gz',
            'WGD_scores': 'WGD_scores.txt.gz',
            'bincov_matrix': f'{self.name}.RD.txt.gz',
            'bincov_matrix_index': f'{self.name}.RD.txt.gz.tbi',
            'bincov_median': 'medianCov.transposed.bed',
            'qc_table': 'evidence_qc_table.tsv',
        }
        for caller in utils.SV_CALLERS:
            for k in ['low', 'high']:
                fname_by_key[f'{caller}_qc_{k}'] = f'{caller}_QC.outlier.{k}'

        return {key: self.get_stage_cohort_prefix(cohort) / fname for key, fname in fname_by_key.items()}

    def queue_jobs(self, cohort: targets.Cohort, inputs: stage.StageInput) -> stage.StageOutput:
        outputs = self.expected_outputs(cohort)

        input_dict = inputs.as_dict_by_target(GatherSampleEvidence)

        jobs = create_evidence_qc_jobs(input_dict=input_dict, output_dict=outputs, cohort=cohort)

        return self.make_outputs(cohort, data=outputs, jobs=jobs)


@stage.stage(
    required_stages=EvidenceQC,
    analysis_type='sv',
)
class CreateSampleBatches(stage.MultiCohortStage):
    """
    uses the values generated in EvidenceQC
    splits the sequencing groups into batches based on median coverage,
    PCR +/- status, and Sex

    The output of this Stage will contain the distinct SG batches to use for the
    following series of Stages. For now, standard practice is to create a separate
    minimal configuration file for each sub-batch, containing the list of SG IDs
    as the `only_sgs` key. The gatk_sv_multisample_1 and gatk_sv_sandwich WFs are
    then run separately for each sub-batch, with the active SGs controlled via the
    config contents.

    The output of this stage is used to generate custom cohorts
    """

    def expected_outputs(self, multicohort: targets.MultiCohort) -> Path:
        return self.prefix / 'sgid_batches.json'

    def queue_jobs(self, multicohort: targets.MultiCohort, inputs: stage.StageInput) -> stage.StageOutput:
        """
        this stage has been clarified - this will only run on Genomes
        Exomes were never a supported case, they have a separate pipeline
        """

        output = self.expected_outputs(multicohort)

        job = create_sample_batches(
            qc_tables=[inputs.as_str(cohort, EvidenceQC, 'qc_table') for cohort in multicohort.get_cohorts()],
            tmp_prefix=self.tmp_prefix,
            output_json=str(output),
        )

        return self.make_outputs(multicohort, data=output, jobs=job)


def cli_main():
    """
    CLI entrypoint - starts up the workflow
    """
    parser = argparse.ArgumentParser()
    parser.add_argument('--dry_run', action='store_true', help='Dry run')
    args = parser.parse_args()

    workflow.run_workflow(name='gatk_sv_singlesample', stages=[CreateSampleBatches], dry_run=args.dry_run)


if __name__ == '__main__':
    cli_main()
