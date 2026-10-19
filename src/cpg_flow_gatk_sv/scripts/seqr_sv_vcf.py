"""
Rewrite a single-dataset subset of the GATK-SV callset into a VCF which seqr's SV loading pipeline accepts

seqr validates the SV callset against a fixed schema, and identifies samples by participant external ID. This:
    - drops excluded samples, and fails if any other expected sample is absent from the VCF
    - replaces SG IDs with participant external IDs
    - adds the header definitions seqr requires: GNOMAD_V4.1_TRUTH_VID, SEQR_INTERNAL_TRUTH_VID and CONC_ST
    - copies the gnomAD v4.1 SV ID from the SVAnnotate external AF field into GNOMAD_V4.1_TRUTH_VID
    - drops records on contigs seqr doesn't load
    - fails on any ALT allele, or CPX_TYPE, which seqr would reject
    - makes variant IDs unique within the file

SEQR_INTERNAL_TRUTH_VID is declared but left empty, which is correct for a first load into seqr.
CONC_ST is declared if absent, values are kept where UpdateStructuralVariantIDs populated them.

The input can be large, so it's streamed line by line. The output is gzip-compressed, not bgzipped.
"""

import gzip
import re
from argparse import ArgumentParser
from collections import Counter, defaultdict
from dataclasses import dataclass, field

import loguru

from cpg_utils.config import config_retrieve

GNOMAD_TRUTH_VID = 'GNOMAD_V4.1_TRUTH_VID'
SEQR_TRUTH_VID = 'SEQR_INTERNAL_TRUTH_VID'

# header definitions seqr requires, which GATK-SV output may lack. Copied from seqr's own test fixture
REQUIRED_HEADER_LINES: dict[tuple[str, str], str] = {
    ('INFO', GNOMAD_TRUTH_VID): (
        f'##INFO=<ID={GNOMAD_TRUTH_VID},Number=1,Type=String,'
        'Description="ID of the matching variant in the gnomAD v4.1 SV sites VCF">'
    ),
    ('INFO', SEQR_TRUTH_VID): (
        f'##INFO=<ID={SEQR_TRUTH_VID},Number=1,Type=String,'
        'Description="ID of the matching variant in the SV callset previously loaded into seqr">'
    ),
    ('FORMAT', 'CONC_ST'): (
        '##FORMAT=<ID=CONC_ST,Number=.,Type=String,Description="The genotype concordance contingency state">'
    ),
}

# seqr skips records on any other contig
SEQR_CONTIGS = frozenset([f'chr{contig}' for contig in [*range(1, 23), 'X', 'Y', 'M']])

# from seqr's loading_pipeline/lib/annotations/enums.py, any other value fails the load
SV_TYPES = frozenset(['BND', 'CPX', 'CTX', 'DEL', 'DUP', 'INS', 'INV', 'CNV'])
SV_TYPE_DETAILS = frozenset(
    [
        'INS_iDEL',
        'INVdel',
        'INVdup',
        'ME',
        'ME:ALU',
        'ME:LINE1',
        'ME:SVA',
        'dDUP',
        'dDUP_iDEL',
        'delINV',
        'delINVdel',
        'delINVdup',
        'dupINV',
        'dupINVdel',
        'dupINVdup',
    ]
)

HEADER_ID_REGEX = re.compile(r'^##(INFO|FORMAT)=<ID=([^,>]+)')
SYMBOLIC_ALT_REGEX = re.compile(r'^<([^<>]+)>$')

# number of fixed VCF columns, CHROM to FORMAT
VCF_FIXED_COLUMNS = 9


@dataclass(frozen=True)
class SampleResolution:
    """The VCF columns to keep, and the participant ID to write for each, in VCF order."""

    column_indices: tuple[int, ...]
    participant_ids: tuple[str, ...]
    families_missing_members: dict[str, list[str]]


@dataclass
class RecordStats:
    total: int = 0
    written: int = 0
    skipped_contigs: Counter = field(default_factory=Counter)
    with_gnomad_vid: int = 0
    renamed_ids: int = 0


def read_sample_map(path: str) -> dict[str, tuple[str, str]]:
    """
    Read the sample map written when the stage was queued

    Args:
        path (str): TSV with no header, columns are SG ID, participant external ID and family ID

    Returns:
        dict of {SG ID: (participant external ID, family ID)}
    """
    sample_map: dict[str, tuple[str, str]] = {}
    with open(path) as handle:
        for line in handle:
            if not line.strip():
                continue
            sg_id, participant_id, family_id = line.rstrip('\n').split('\t')
            sample_map[sg_id] = (participant_id, family_id)
    return sample_map


def read_exclusions(path: str) -> set[str]:
    """Read the combined exclusion list, one SG ID per line."""
    with open(path) as handle:
        return {line.strip() for line in handle if line.strip()}


def resolve_samples(
    vcf_samples: list[str],
    sample_map: dict[str, tuple[str, str]],
    exclusions: set[str],
) -> SampleResolution:
    """
    Work out which VCF columns to keep, and what to call each one

    Excluded SGs are dropped, as the pipeline already removed them from consideration
    Any other SG in the sample map which is absent from the VCF is an error, as seqr can't load it
    Any SG in the VCF which isn't in the sample map is dropped, it belongs to another dataset or family

    Args:
        vcf_samples (list[str]): sample IDs from the #CHROM line, in order
        sample_map (dict): {SG ID: (participant external ID, family ID)} for every SG to export
        exclusions (set[str]): SG IDs removed by the pipeline's QC

    Returns:
        SampleResolution
    """
    expected = set(sample_map) - exclusions
    if missing := sorted(expected - set(vcf_samples)):
        raise ValueError(f'{len(missing)} SG IDs expected in the VCF are absent: {", ".join(missing)}')

    kept = [(index, sg_id) for index, sg_id in enumerate(vcf_samples) if sg_id in expected]
    if not kept:
        raise ValueError('No samples left to export after removing exclusions')

    participant_ids = tuple(sample_map[sg_id][0] for _index, sg_id in kept)
    if duplicates := sorted(pid for pid, count in Counter(participant_ids).items() if count > 1):
        raise ValueError(f'Participant external IDs shared by multiple SGs: {", ".join(duplicates)}')

    # the rest of the family still loads, but seqr will refuse the load if an excluded member already has SV data
    families_missing_members: dict[str, list[str]] = defaultdict(list)
    for sg_id in sorted(set(sample_map) & exclusions):
        participant_id, family_id = sample_map[sg_id]
        families_missing_members[family_id].append(f'{participant_id} ({sg_id})')

    return SampleResolution(
        column_indices=tuple(index for index, _sg_id in kept),
        participant_ids=participant_ids,
        families_missing_members=dict(families_missing_members),
    )


def missing_header_lines(header_lines: list[str]) -> list[str]:
    """Return the required header definitions which aren't already in the header."""
    present: set[tuple[str, str]] = set()
    for line in header_lines:
        if match := HEADER_ID_REGEX.match(line):
            present.add((match.group(1), match.group(2)))
    return [definition for key, definition in REQUIRED_HEADER_LINES.items() if key not in present]


def info_value(info: str, key: str) -> str | None:
    """Get the value of a key=value INFO entry, None if it's absent."""
    prefix = f'{key}='
    for entry in info.split(';'):
        if entry.startswith(prefix):
            return entry[len(prefix) :]
    return None


def validate_alleles(chrom: str, pos: str, variant_id: str, alt: str, info: str) -> None:
    """
    Fail on anything seqr would reject: non-symbolic ALTs, unknown SV types, unknown INS subtypes or CPX_TYPEs
    seqr only checks the subtype on INS alleles, e.g. <INS:ME:ALU>, and checks CPX_TYPE for CPX
    """
    location = f'{variant_id} at {chrom}:{pos}'
    for allele in alt.split(','):
        if not (match := SYMBOLIC_ALT_REGEX.match(allele)):
            raise ValueError(f'{location} has a non-symbolic ALT allele: {allele}')
        sv_type, _, detail = match.group(1).partition(':')
        if sv_type not in SV_TYPES:
            raise ValueError(f'{location} has an SV type seqr does not support: {allele}')
        if sv_type == 'INS' and detail and detail not in SV_TYPE_DETAILS:
            raise ValueError(f'{location} has an INS subtype seqr does not support: {allele}')
        if sv_type == 'CPX' and (cpx_type := info_value(info, 'CPX_TYPE')) and cpx_type not in SV_TYPE_DETAILS:
            raise ValueError(f'{location} has a CPX_TYPE seqr does not support: {cpx_type}')


def unique_id(variant_id: str, seen_ids: set[str]) -> str:
    """
    Return the ID, or if it's already used in this file, the ID with the lowest free numeric suffix
    e.g. DEL_1-100-200, DEL_1-100-200_2, DEL_1-100-200_3
    """
    if variant_id not in seen_ids:
        return variant_id
    suffix = 2
    while f'{variant_id}_{suffix}' in seen_ids:
        suffix += 1
    return f'{variant_id}_{suffix}'


def transform_record(
    line: str,
    column_indices: tuple[int, ...],
    gnomad_source_field: str,
    seen_ids: set[str],
    stats: RecordStats,
) -> str | None:
    """
    Rewrite a single VCF record for seqr, or return None if seqr won't load it

    Adds the record's (possibly renamed) ID to seen_ids, and updates stats
    """
    stats.total += 1
    columns = line.rstrip('\n').split('\t')
    chrom, pos, variant_id, _ref, alt, _qual, _filters, info = columns[:8]

    if chrom not in SEQR_CONTIGS:
        stats.skipped_contigs[chrom] += 1
        return None

    if variant_id in {'', '.'}:
        raise ValueError(f'Record at {chrom}:{pos} has no variant ID, which seqr requires')

    # run validation - if we encounter failures in a permissive mode, log instead of failing
    try:
        validate_alleles(chrom, pos, variant_id, alt, info)
    except ValueError as ve:
        if config_retrieve(['workflow', 'sv_transform_permissive'], False):
            loguru.logger.warning(ve)
            return None
        raise ve

    new_id = unique_id(variant_id, seen_ids)
    seen_ids.add(new_id)
    if new_id != variant_id:
        stats.renamed_ids += 1

    if info_value(info, GNOMAD_TRUTH_VID) is None and (gnomad_vid := info_value(info, gnomad_source_field)):
        entry = f'{GNOMAD_TRUTH_VID}={gnomad_vid}'
        info = entry if info in {'', '.'} else f'{info};{entry}'
    if info_value(info, GNOMAD_TRUTH_VID):
        stats.with_gnomad_vid += 1

    samples = [columns[VCF_FIXED_COLUMNS + index] for index in column_indices]
    stats.written += 1
    return '\t'.join([chrom, pos, new_id, *columns[3:7], info, columns[8], *samples]) + '\n'


def log_summary(stats: RecordStats, resolution: SampleResolution) -> None:
    loguru.logger.info(f'Wrote {stats.written}/{stats.total} records for {len(resolution.participant_ids)} samples')
    for chrom, count in sorted(stats.skipped_contigs.items()):
        loguru.logger.info(f'Skipped {count} records on {chrom}, seqr does not load this contig')
    if stats.renamed_ids:
        loguru.logger.warning(f'{stats.renamed_ids} records had an ID already used in this file, a suffix was added')
    loguru.logger.info(f'{stats.with_gnomad_vid}/{stats.written} records have a {GNOMAD_TRUTH_VID}')
    for family_id, members in sorted(resolution.families_missing_members.items()):
        loguru.logger.warning(
            f'Family {family_id} has members excluded by QC, who will have no SV calls: {", ".join(members)}'
        )


def main(input_vcf: str, sample_map_file: str, exclusions_file: str, gnomad_source_field: str, output: str) -> None:
    """
    Args:
        input_vcf (str): gzipped VCF, already subset to (roughly) this dataset's samples
        sample_map_file (str): TSV of SG ID, participant external ID, family ID, for every SG to export
        exclusions_file (str): SG IDs removed by the pipeline's QC, one per line
        gnomad_source_field (str): INFO field holding the gnomAD v4.1 SV ID, written by SVAnnotate
        output (str): path to write the gzipped VCF to
    """
    sample_map = read_sample_map(sample_map_file)
    exclusions = read_exclusions(exclusions_file)

    stats = RecordStats()
    seen_ids: set[str] = set()
    header_lines: list[str] = []
    resolution: SampleResolution | None = None

    # fast compression, bcftools recompresses this as bgzip
    with gzip.open(input_vcf, 'rt') as f_in, gzip.open(output, 'wt', compresslevel=1) as f_out:
        for line in f_in:
            if line.startswith('##'):
                header_lines.append(line.rstrip('\n'))
                continue

            if line.startswith('#CHROM'):
                columns = line.rstrip('\n').split('\t')
                resolution = resolve_samples(columns[VCF_FIXED_COLUMNS:], sample_map, exclusions)
                header_lines.extend(missing_header_lines(header_lines))
                header_lines.append('\t'.join([*columns[:VCF_FIXED_COLUMNS], *resolution.participant_ids]))
                f_out.write('\n'.join(header_lines) + '\n')
                continue

            if resolution is None:
                raise ValueError(f'{input_vcf} has a record before the #CHROM header line')

            if new_line := transform_record(line, resolution.column_indices, gnomad_source_field, seen_ids, stats):
                f_out.write(new_line)

    if resolution is None:
        raise ValueError(f'{input_vcf} has no #CHROM header line')
    if not stats.written:
        raise ValueError(f'No records written from {input_vcf}, check it uses chr-prefixed GRCh38 contigs')

    log_summary(stats, resolution)


def cli_main():
    parser = ArgumentParser(description='Rewrite a single-dataset SV VCF into a format seqr can load')
    parser.add_argument('--input', required=True, help='gzipped VCF, subset to this dataset')
    parser.add_argument('--sample_map', required=True, help='TSV of SG ID, participant external ID, family ID')
    parser.add_argument('--exclusions', required=True, help='SG IDs removed by QC, one per line')
    parser.add_argument('--gnomad_field', required=True, help='INFO field holding the gnomAD v4.1 SV ID')
    parser.add_argument('--output', required=True, help='path to write the gzipped VCF to')
    args = parser.parse_args()
    main(
        input_vcf=args.input,
        sample_map_file=args.sample_map,
        exclusions_file=args.exclusions,
        gnomad_source_field=args.gnomad_field,
        output=args.output,
    )


if __name__ == '__main__':
    cli_main()
