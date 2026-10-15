import gzip

from cpg_flow_gatk_sv.scripts.rename_sv_ids import make_coordinate_id, rename_sv_ids

HEADER = '##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n'


def record(pos: str, end: str, truth_vid: str | None = None, alt: str = '<DEL>') -> str:
    info = f'END={end};SVTYPE=DEL'
    if truth_vid:
        info = f'{info};TRUTH_VID={truth_vid}'
    return '\t'.join(['chr1', pos, 'all_samples.chr1.final_cleanup_DEL_1', 'N', alt, '999', 'PASS', info]) + '\n'


def run(tmp_path, records: list[str], skip_prior_names: bool = False) -> list[str]:
    input_vcf = tmp_path / 'input.vcf.gz'
    output_vcf = tmp_path / 'output.vcf.gz'
    with gzip.open(input_vcf, 'wt') as handle:
        handle.write(HEADER)
        handle.writelines(records)
    rename_sv_ids(str(input_vcf), str(output_vcf), skip_prior_names=skip_prior_names)
    with gzip.open(output_vcf, 'rt') as handle:
        return [line.split('\t')[2] for line in handle if not line.startswith('#')]


def test_previous_id_is_used(tmp_path):
    assert run(tmp_path, [record('1005', '2003', truth_vid='DEL_1-1000-2000')]) == ['DEL_1-1000-2000']


def test_previous_id_is_ignored_with_skip_prior_names(tmp_path):
    ids = run(tmp_path, [record('1005', '2003', truth_vid='DEL_1-1000-2000')], skip_prior_names=True)
    assert ids == ['DEL_1-1005-2003']


def test_shared_previous_id_is_only_taken_once(tmp_path):
    ids = run(
        tmp_path,
        [
            record('1005', '2003', truth_vid='DEL_1-1000-2000'),
            record('1010', '1995', truth_vid='DEL_1-1000-2000'),
            record('5000', '6000'),
        ],
    )
    assert ids == ['DEL_1-1000-2000', 'DEL_1-1010-1995', 'DEL_1-5000-6000']


def test_previous_id_taken_by_earlier_coordinate_id_is_not_reused(tmp_path):
    ids = run(tmp_path, [record('1000', '2000'), record('1005', '2003', truth_vid='DEL_1-1000-2000')])
    assert ids == ['DEL_1-1000-2000', 'DEL_1-1005-2003']


def test_make_coordinate_id_formats():
    bnd = ['chr1', '100', '.', 'N', '<BND>', '.', 'PASS', '']
    assert make_coordinate_id(bnd, {'CHR2': 'chr5', 'END2': '200', 'END': '100'}) == 'BND_1-100_5-200'
    ins = ['chr1', '100', '.', 'N', '<INS:ME:ALU>', '.', 'PASS', '']
    assert make_coordinate_id(ins, {'SVLEN': '300', 'END': '101'}) == 'INS:ME:ALU_1-100_ins300'
    dup = ['chrX', '100', '.', 'N', '<DUP>', '.', 'PASS', '']
    assert make_coordinate_id(dup, {'END': '500'}) == 'DUP_X-100-500'


def test_exact_coordinate_match_takes_previous_id_over_earlier_match(tmp_path):
    ids = run(
        tmp_path,
        [
            record('990', '2000', truth_vid='DEL_1-1000-2000'),
            record('1000', '2000', truth_vid='DEL_1-1000-2000'),
        ],
    )
    assert ids == ['DEL_1-990-2000', 'DEL_1-1000-2000']
