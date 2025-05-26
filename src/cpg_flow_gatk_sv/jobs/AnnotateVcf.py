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
        'external_af_population': config.config_retrieve(['references', 'gatk_sv', 'external_af_population']),
        'external_af_ref_prefix': config.config_retrieve(['references', 'gatk_sv', 'external_af_ref_bed_prefix']),
        'external_af_ref_bed': config.config_retrieve(['references', 'gnomad_sv']),
        'use_hail': False,
    }

    input_dict |= utils.get_references(
        [
            'noncoding_bed',
            'protein_coding_gtf',
            {'contig_list': 'primary_contigs_list'},
        ],
    )

    input_dict |= {
        key: config.config_retrieve(['images', key])
        for key in [
            'sv_pipeline_docker',
            'sv_base_mini_docker',
            'gatk_docker',
        ]
    }
    return utils.add_gatk_sv_jobs(
        dataset=multicohort.analysis_dataset,
        wfl_name='AnnotateVcf',
        input_dict=input_dict,
        expected_out_dict=outputs,
        labels={'stage': 'annotatevcf', config.AR_GUID_NAME: config.try_get_ar_guid()},
    )
