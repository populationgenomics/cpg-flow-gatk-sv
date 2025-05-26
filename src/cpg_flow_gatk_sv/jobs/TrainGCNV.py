import random
from typing import Any

from cpg_flow import targets
from cpg_flow_gatk_sv import utils
from cpg_utils import Path, config, to_path


def add_train_gcnv_jobs(cohort: targets.Cohort, output_dict: dict[str, Path]):
    """

    Returns:

    """

    # optionally do sample subsetting
    sample_n = config.config_retrieve(['train_gcnv', 'sample_n'], 100)
    all_sgids = cohort.get_sequencing_groups()

    # TODO force a minimum number of samples?
    sgs_sampled_from_cohort = random.sample(all_sgids, min(len(all_sgids), sample_n))

    # new fasta getter method
    fasta_file = utils.get_fasta_string()

    # pull in the basic input dict
    cromwell_input_dict: dict[str, Any] = {
        'cohort': cohort.id,
        'samples': [sg.id for sg in sgs_sampled_from_cohort],
        'count_files': [
            str(sequencing_group.make_sv_evidence_path / f'{sequencing_group.id}.coverage_counts.tsv.gz')
            for sequencing_group in sgs_sampled_from_cohort
        ],
        'ref_copy_number_autosomal_contigs': 2,
        'allosomal_contigs': ['chrX', 'chrY'],
        'reference_fasta': fasta_file,
        'reference_index': f'{fasta_file}.fai',
        'reference_dict': to_path(fasta_file).with_suffix('.dict'),
        'contig_ploidy_priors': config.reference_path('gatk_sv/contig_ploidy_priors'),
        'num_intervals_per_scatter': 5000,
    }

    # add the images required for this step
    cromwell_input_dict |= utils.get_images(
        [
            'sv_base_mini_docker',
            'linux_docker',
            'gatk_docker',
            'condense_counts_docker',
            'sv_pipeline_docker',
        ]
    )

    # billing labels!
    # https://cromwell.readthedocs.io/en/stable/wf_options/Google/
    # these must conform to the regex [a-z]([-a-z0-9]*[a-z0-9])?
    billing_labels = {
        'stage': 'traingcnv',
        config.AR_GUID_NAME: config.try_get_ar_guid(),
    }

    return utils.add_gatk_sv_jobs(
        dataset=cohort.dataset,
        wfl_name='TrainGCNV',
        input_dict=cromwell_input_dict,
        expected_out_dict=output_dict,
        labels=billing_labels,
        job_size=utils.CromwellJobSizes.MEDIUM,
    )
