"""
method(s) for GATK SV PythonJobs

The ~1GB VCF file buffs up to ~55GB when uncompressed
We need to be careful to read and write compressed versions
And don't store any substantial amount of data in memory
"""

from argparse import ArgumentParser


def rename_sv_ids(input_vcf: str, output_vcf: str, skip_prior_names: bool = False):
    """
    A Python method to call as a PythonJob, edits content of the VCF
    Replaces the standard ID with a guaranteed-unique ID (hopefully)
    The new ID compounds Type, Chr & Start, then either end or Chr2 & End2

    Writes the file back out to the specified path

    Args:
        input_vcf (str): path to temp file generated by merging
        output_vcf (str): path to write uncompressed edited version to
        skip_prior_names (bool): if True, don't use the TRUTH_VID ID
    """
    import gzip

    # crack open that VCF and have a little sip
    with gzip.open(input_vcf, 'rt') as f, gzip.open(output_vcf, 'wt') as f_out:
        for line in f:
            # don't alter current header lines
            if line.startswith('#'):
                f_out.write(line)
                continue

            # split on tabs
            l_split = line.split('\t')

            # some entries are not key-value, so skip them, e.g.
            # AN_Orig=61;END=56855888;SVTYPE=DUP;BOTHSIDES_SUPPORT;TRUTH_VID=CNV_1-46855888-56855888  # noqa: ERA001
            info_dict: dict[str, str] = {}
            for entry in l_split[7].split(';'):
                if '=' in entry:
                    key, value = entry.split('=')
                    info_dict[key] = value

            # if this matches a previous ID, use that. TRUTH_VID = Truth dataset Variant ID
            # In this context 'truth' was previous results for consistency. Not a validation truth set
            # boolean here indicates when we should ignore TRUTH_VID
            # this occurs when we have not yet created meaningful IDs for this cohort
            if info_dict.get('TRUTH_VID') and not skip_prior_names:
                l_split[2] = info_dict['TRUTH_VID']
                f_out.write('\t'.join(l_split))
                continue

            # otherwise generate a new ID
            # strip out the "chr" prefix to abbreviate String
            chrom = l_split[0].removeprefix('chr')
            start = l_split[1]

            # e.g. <DEL> -> DEL
            alt_allele = l_split[4][1:-1]

            # maybe there's a second chromosome?
            if all(key in info_dict for key in ['CHR2', 'END2']):
                chrom2 = info_dict['CHR2'].removeprefix('chr')
                end2 = info_dict['END2']
                # e.g. BND_1-12345_2-67890
                key = f'{alt_allele}_{chrom}-{start}_{chrom2}-{end2}'

            elif 'INS' in alt_allele and 'SVLEN' in info_dict:
                # e.g. INS:ME:ALU_1-12345-67890
                key = f'{alt_allele}_{chrom}-{start}_ins{info_dict["SVLEN"]}'

            else:
                # e.g. CNV_1-12345-67890
                key = f'{alt_allele}_{chrom}-{start}-{info_dict["END"]}'

            # some manta calls are exactly identical...
            # we can make them unique by shoving stuff on the end, but we'll have issues
            # matching exactly to this call in future
            l_split[2] = key

            # rebuild the line
            f_out.write('\t'.join(l_split))


if __name__ == '__main__':
    parser = ArgumentParser(description='Rename SV IDs in a VCF')
    parser.add_argument('--input_vcf', help='path to temp file generated by merging')
    parser.add_argument('--output_vcf', help='path to write uncompressed edited version to')
    parser.add_argument(
        '--skip_prior_names',
        action='store_true',
        help='if True, do not use the TRUTH_VID ID from the VCF',
    )
    args = parser.parse_args()

    rename_sv_ids(args.input_vcf, args.output_vcf, skip_prior_names=args.skip_prior_names)
