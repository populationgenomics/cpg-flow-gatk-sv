from unittest.mock import MagicMock, patch

import hailtop.batch as hb

from cpg_flow_gatk_sv.jobs import SeqrSvVcf

SG_TO_PARTICIPANT = {'SG1': 'P1', 'SG2': 'P2', 'SG3': 'P3'}


def test_split_job_writes_each_exported_participant_separately(tmp_path):
    batch = hb.Batch(backend=hb.LocalBackend())
    split_outputs = {'SG1': tmp_path / 'split' / 'P1.vcf.gz', 'SG3': tmp_path / 'split' / 'P3.vcf.gz'}
    with (
        patch.object(SeqrSvVcf.hail_batch, 'get_batch', return_value=batch),
        patch.object(SeqrSvVcf.config, 'config_retrieve', return_value='image'),
        patch.object(SeqrSvVcf, 'write_sample_files'),
    ):
        jobs = SeqrSvVcf.create_seqr_sv_vcf_jobs(
            dataset=MagicMock(),
            input_vcf=str(tmp_path / 'callset.vcf.bgz'),
            sg_to_participant=SG_TO_PARTICIPANT,
            sample_map_file=tmp_path / 'sample_map.tsv',
            sg_id_file=tmp_path / 'sg_ids.txt',
            exclusions=str(tmp_path / 'exclusions.txt'),
            output_vcf=tmp_path / 'seqr.vcf.bgz',
            split_outputs=split_outputs,
            job_attrs={},
        )

    split_job = next(job for job in jobs if 'split by participant' in job.name)
    commands = '\n'.join(split_job._command)  # noqa: SLF001
    assert '"split/P1.vcf.gz"' in commands
    assert '"split/P3.vcf.gz"' in commands
    assert 'P2' not in commands
    for sg_id, output in split_outputs.items():
        assert split_job[sg_id]['vcf.gz']._output_paths == {str(output)}  # noqa: SLF001
        assert split_job[sg_id]['vcf.gz.tbi']._output_paths == {f'{output}.tbi'}  # noqa: SLF001
