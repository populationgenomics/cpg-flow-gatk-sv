from typing import TYPE_CHECKING

from cpg_flow import targets
from cpg_flow_gatk_sv import utils
from cpg_utils import Path, config

if TYPE_CHECKING:
    from hailtop.batch.job import BashJob


def create_svannotate_jobs(
    multicohort: targets.MultiCohort,
    input_vcf: str,
    pedigree: str,
    outputs: dict[str, Path],
) -> list['BashJob']:
    input_dict = {
        'vcf': input_vcf,
        'prefix': multicohort.name,
        'ped_file': pedigree,
        'sv_per_shard': 5000,
        'external_af_population': config.config_retrieve(['references', 'external_af_population']),
        'external_af_ref_prefix': config.config_retrieve(['references', 'external_af_ref_bed_prefix']),
        'external_af_ref_bed': config.config_retrieve(['references', 'gnomad_sv']),
        'use_hail': False,
        'noncoding_bed': config.config_retrieve(['references', 'noncoding_bed']),
        'protein_coding_gtf': config.config_retrieve(['references', 'protein_coding_gtf']),
        'contig_list': config.config_retrieve(['references', 'primary_contigs_list']),
        'gatk_docker': config.config_retrieve(['images', 'gq_recalibrator_docker']),
        'sv_base_mini_docker': config.config_retrieve(['images', 'sv_base_mini_docker']),
        'sv_pipeline_docker': config.config_retrieve(['images', 'sv_pipeline_docker']),
    }

    return utils.add_gatk_sv_jobs(
        dataset=multicohort.analysis_dataset,
        wfl_name='AnnotateVcf',
        input_dict=input_dict,
        expected_out_dict=outputs,
        labels={'stage': 'annotatevcf', config.AR_GUID_NAME: config.try_get_ar_guid()},
    )
