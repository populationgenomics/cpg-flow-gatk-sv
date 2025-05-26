from typing import TYPE_CHECKING

from cpg_flow import targets
from cpg_flow_gatk_sv import utils
from cpg_utils import Path, config, to_path

if TYPE_CHECKING:
    from hailtop.batch.job import BashJob


def create_sv_concordance_jobs(
    multicohort: targets.MultiCohort,
    formatvcf_output: str,
    joinrawcalls_output: str,
    outputs: dict[str, Path],
) -> list['BashJob']:
    fasta_file = utils.get_fasta_string()

    input_dict = {
        'output_prefix': multicohort.name,
        'reference_dict': str(to_path(fasta_file).with_suffix('.dict')),
        'eval_vcf': formatvcf_output,
        'truth_vcf': joinrawcalls_output,
    }

    # add the images required for this step
    input_dict |= {
        key: config.config_retrieve(['images', key])
        for key in [
            'gatk_docker',
            'sv_base_mini_docker',
        ]
    }
    input_dict |= utils.get_references(
        [
            {'contig_list': 'primary_contigs_list'},
        ]
    )

    return utils.add_gatk_sv_jobs(
        dataset=multicohort.analysis_dataset,
        wfl_name='SVConcordance',
        input_dict=input_dict,
        expected_out_dict=outputs,
        labels={'stage': 'svconcordance', config.AR_GUID_NAME: config.try_get_ar_guid()},
    )
