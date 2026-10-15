import gzip

import pytest

from cpg_flow_gatk_sv.scripts.seqr_sv_vcf import (
    GNOMAD_TRUTH_VID,
    REQUIRED_HEADER_LINES,
    RecordStats,
    main,
    missing_header_lines,
    resolve_samples,
    transform_record,
    unique_id,
    validate_alleles,
)

GNOMAD_FIELD = 'gnomad_v4.1_sv_SVID'
FORMAT = 'GT:GQ:RD_CN'

SAMPLE_MAP = {
    'SG1': ('P1', 'FAM1'),
    'SG2': ('P2', 'FAM1'),
    'SG3': ('P3', 'FAM2'),
}

HEADER = [
    '##fileformat=VCFv4.2',
    '##contig=<ID=chr1,length=248956422>',
    '##INFO=<ID=END,Number=1,Type=Integer,Description="End position">',
    f'##INFO=<ID={GNOMAD_FIELD},Number=1,Type=String,Description="gnomAD SV ID">',
    '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
    '##FORMAT=<ID=CONC_ST,Number=.,Type=String,Description="The genotype concordance contingency state">',
]


def record(chrom='chr1', pos='100', vid='DEL_1-100-200', alt='<DEL>', info='END=200', samples=('0/1:99:1',)) -> str:
    return '\t'.join([chrom, pos, vid, 'N', alt, '999', 'PASS', info, FORMAT, *samples]) + '\n'


def write_vcf(path, samples: list[str], records: list[str]) -> None:
    with gzip.open(path, 'wt') as handle:
        handle.write('\n'.join(HEADER) + '\n')
        handle.write('\t'.join(['#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT', *samples]))
        handle.write('\n')
        handle.writelines(records)


def test_resolve_samples_renames_in_vcf_order_and_drops_unmapped():
    resolution = resolve_samples(['SG3', 'OTHER', 'SG1', 'SG2'], SAMPLE_MAP, exclusions=set())
    assert resolution.column_indices == (0, 2, 3)
    assert resolution.participant_ids == ('P3', 'P1', 'P2')
    assert resolution.families_missing_members == {}


def test_resolve_samples_drops_exclusions_and_reports_their_families():
    resolution = resolve_samples(['SG1', 'SG3'], SAMPLE_MAP, exclusions={'SG2'})
    assert resolution.participant_ids == ('P1', 'P3')
    assert resolution.families_missing_members == {'FAM1': ['P2 (SG2)']}


def test_resolve_samples_drops_excluded_sample_present_in_vcf():
    resolution = resolve_samples(['SG1', 'SG2', 'SG3'], SAMPLE_MAP, exclusions={'SG2'})
    assert resolution.column_indices == (0, 2)


def test_resolve_samples_fails_on_missing_sample():
    with pytest.raises(ValueError, match='absent: SG2'):
        resolve_samples(['SG1', 'SG3'], SAMPLE_MAP, exclusions=set())


def test_resolve_samples_fails_on_shared_participant():
    sample_map = SAMPLE_MAP | {'SG4': ('P1', 'FAM1')}
    with pytest.raises(ValueError, match='shared by multiple SGs: P1'):
        resolve_samples(['SG1', 'SG2', 'SG3', 'SG4'], sample_map, exclusions=set())


def test_resolve_samples_fails_when_nothing_left():
    with pytest.raises(ValueError, match='No samples left'):
        resolve_samples(['SG1'], {'SG1': ('P1', 'FAM1')}, exclusions={'SG1'})


def test_missing_header_lines_only_adds_absent_definitions():
    added = missing_header_lines(HEADER)
    assert added == [
        REQUIRED_HEADER_LINES[('INFO', GNOMAD_TRUTH_VID)],
        REQUIRED_HEADER_LINES[('INFO', 'SEQR_INTERNAL_TRUTH_VID')],
    ]


@pytest.mark.parametrize(
    ('alt', 'info'),
    [
        ('<DEL>', 'END=200'),
        ('<INS:ME:ALU>', 'END=200'),
        ('<INS:ME>', 'END=200'),
        ('<CPX>', 'CPX_TYPE=dupINVdel'),
        ('<BND>', 'CHR2=chr5'),
        ('<CNV>', 'END=200'),
        ('<CTX>', 'END=200'),
    ],
)
def test_validate_alleles_accepts_seqr_types(alt, info):
    validate_alleles('chr1', '100', 'ID', alt, info)


@pytest.mark.parametrize(
    ('alt', 'info', 'message'),
    [
        ('A', 'END=200', 'non-symbolic'),
        ('<DEL>,T', 'END=200', 'non-symbolic'),
        ('<CN0>', 'END=200', 'SV type'),
        ('<INS:UNK>', 'END=200', 'INS subtype'),
        ('<CPX>', 'CPX_TYPE=madeUp', 'CPX_TYPE'),
    ],
)
def test_validate_alleles_rejects_unsupported(alt, info, message):
    with pytest.raises(ValueError, match=message):
        validate_alleles('chr1', '100', 'ID', alt, info)


def test_unique_id_takes_lowest_free_suffix():
    assert unique_id('A', set()) == 'A'
    assert unique_id('A', {'A'}) == 'A_2'
    assert unique_id('A', {'A', 'A_2', 'A_3'}) == 'A_4'


def test_transform_record_copies_gnomad_id_and_subsets_samples():
    stats = RecordStats()
    line = record(info=f'END=200;{GNOMAD_FIELD}=gnomAD-SV_v3_DEL_chr1_abc', samples=('0/1:99:1', '0/0:99:2'))
    result = transform_record(line, (1,), GNOMAD_FIELD, set(), stats)
    columns = result.rstrip('\n').split('\t')
    assert (
        columns[7] == f'END=200;{GNOMAD_FIELD}=gnomAD-SV_v3_DEL_chr1_abc;{GNOMAD_TRUTH_VID}=gnomAD-SV_v3_DEL_chr1_abc'
    )
    assert columns[9:] == ['0/0:99:2']
    assert stats.with_gnomad_vid == 1


def test_transform_record_keeps_existing_gnomad_truth_vid():
    stats = RecordStats()
    info = f'END=200;{GNOMAD_FIELD}=other;{GNOMAD_TRUTH_VID}=kept'
    result = transform_record(record(info=info), (0,), GNOMAD_FIELD, set(), stats)
    assert result.split('\t')[7] == info


def test_transform_record_skips_other_contigs():
    stats = RecordStats()
    assert transform_record(record(chrom='chrUn_KI270742v1'), (0,), GNOMAD_FIELD, set(), stats) is None
    assert transform_record(record(chrom='1'), (0,), GNOMAD_FIELD, set(), stats) is None
    assert stats.skipped_contigs == {'chrUn_KI270742v1': 1, '1': 1}
    assert stats.written == 0


def test_transform_record_fails_without_id():
    with pytest.raises(ValueError, match='no variant ID'):
        transform_record(record(vid='.'), (0,), GNOMAD_FIELD, set(), RecordStats())


def test_main_writes_seqr_ready_vcf(tmp_path):
    input_vcf = tmp_path / 'input.vcf.gz'
    write_vcf(
        input_vcf,
        samples=['SG1', 'SG3'],
        records=[
            record(info=f'END=200;{GNOMAD_FIELD}=gnomAD-SV_v3_DEL_chr1_abc', samples=('0/1:99:1', '0/0:99:2')),
            record(samples=('0/0:99:2', '0/1:99:1')),
            record(chrom='chrUn_KI270742v1', vid='DEL_Un-1-2', samples=('0/0:99:2', '0/0:99:2')),
            record(pos='300', vid='DUP_1-300-400', alt='<DUP>', info='END=400', samples=('0/0:99:2', '1/1:99:4')),
        ],
    )
    sample_map = tmp_path / 'sample_map.tsv'
    sample_map.write_text(''.join(f'{sg}\t{pid}\t{fam}\n' for sg, (pid, fam) in SAMPLE_MAP.items()))
    exclusions = tmp_path / 'exclusions.txt'
    exclusions.write_text('SG2\n')
    output = tmp_path / 'output.vcf.gz'

    main(str(input_vcf), str(sample_map), str(exclusions), GNOMAD_FIELD, str(output))

    with gzip.open(output, 'rt') as handle:
        lines = handle.read().splitlines()
    header = [line for line in lines if line.startswith('##')]
    chrom_line = next(line for line in lines if line.startswith('#CHROM'))
    records = [line.split('\t') for line in lines if not line.startswith('#')]

    assert all(definition in header for definition in REQUIRED_HEADER_LINES.values())
    assert sum(line.startswith('##FORMAT=<ID=CONC_ST,') for line in header) == 1
    assert chrom_line.split('\t')[9:] == ['P1', 'P3']
    assert [columns[2] for columns in records] == ['DEL_1-100-200', 'DEL_1-100-200_2', 'DUP_1-300-400']
    assert f'{GNOMAD_TRUTH_VID}=gnomAD-SV_v3_DEL_chr1_abc' in records[0][7]
    assert GNOMAD_TRUTH_VID not in records[1][7]


def test_main_fails_when_no_records_on_seqr_contigs(tmp_path):
    input_vcf = tmp_path / 'input.vcf.gz'
    write_vcf(input_vcf, samples=['SG1'], records=[record(chrom='1')])
    sample_map = tmp_path / 'sample_map.tsv'
    sample_map.write_text('SG1\tP1\tFAM1\n')
    exclusions = tmp_path / 'exclusions.txt'
    exclusions.write_text('')

    with pytest.raises(ValueError, match='No records written'):
        main(str(input_vcf), str(sample_map), str(exclusions), GNOMAD_FIELD, str(tmp_path / 'out.vcf.gz'))
