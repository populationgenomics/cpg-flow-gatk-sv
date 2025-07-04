from typing import TYPE_CHECKING

from cpg_flow import targets
from cpg_flow_gatk_sv import utils
from cpg_utils import Path, config, to_path

if TYPE_CHECKING:
    from hailtop.batch.job import BashJob


def create_genotypebatch_jobs(
    cohort: targets.Cohort,
    gatherbatchevidence_outputs: dict[str, Path],
    filterbatch_outputs: dict[str, Path],
    mergebatchsites_outputs: dict[str, Path],
    outputs: dict[str, Path],
) -> list['BashJob']:
    fasta_file = utils.get_fasta_string()

    input_dict = {
        'batch': cohort.id,
        'n_per_split': config.config_retrieve(
            ['resource_overrides', 'GenotypeBatch', 'n_per_split'],
            5000,
        ),
        'n_RD_genotype_bins': config.config_retrieve(
            ['resource_overrides', 'GenotypeBatch', 'n_RD_genotype_bins'],
            100000,
        ),
        'coveragefile': gatherbatchevidence_outputs['merged_bincov'],
        'coveragefile_index': gatherbatchevidence_outputs['merged_bincov_index'],
        'discfile': gatherbatchevidence_outputs['merged_PE'],
        'discfile_index': gatherbatchevidence_outputs['merged_PE_index'],
        'splitfile': gatherbatchevidence_outputs['merged_SR'],
        'splitfile_index': gatherbatchevidence_outputs['merged_SR_index'],
        'medianfile': gatherbatchevidence_outputs['median_cov'],
        'rf_cutoffs': filterbatch_outputs['cutoffs'],
        'ref_dict': str(to_path(fasta_file).with_suffix('.dict')),
        'reference_build': 'hg38',
        'batch_depth_vcf': filterbatch_outputs['filtered_depth_vcf'],
        'batch_pesr_vcf': filterbatch_outputs['filtered_pesr_vcf'],
        'cohort_depth_vcf': mergebatchsites_outputs['cohort_depth_vcf'],
        'cohort_pesr_vcf': mergebatchsites_outputs['cohort_pesr_vcf'],
        'primary_contigs_list': config.config_retrieve(['references', 'primary_contigs_list']),
        'seed_cutoffs': config.config_retrieve(['references', 'seed_cutoffs']),
        'pesr_exclude_list': config.config_retrieve(['references', 'pesr_exclude_list']),
        'bin_exclude': config.config_retrieve(['references', 'bin_exclude']),
        'linux_docker': config.config_retrieve(['images', 'linux_docker']),
        'sv_base_mini_docker': config.config_retrieve(['images', 'sv_base_mini_docker']),
        'sv_pipeline_docker': config.config_retrieve(['images', 'sv_pipeline_docker']),
    }

    return utils.add_gatk_sv_jobs(
        dataset=cohort.dataset,
        wfl_name='GenotypeBatch',
        input_dict=input_dict,
        expected_out_dict=outputs,
        labels={'stage': 'genotypebatch', config.AR_GUID_NAME: config.try_get_ar_guid()},
    )
