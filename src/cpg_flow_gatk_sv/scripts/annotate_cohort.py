"""
Hail Query functions for seqr loader; SV edition.
"""

import argparse
import gzip

import loguru

import hail as hl

from cpg_utils import config, hail_batch

# I'm just going to go ahead and steal these constants from their seqr loader
BOTHSIDES_SUPPORT = 'BOTHSIDES_SUPPORT'
GENE_SYMBOL = 'gene_symbol'
GENE_ID = 'gene_id'
MAJOR_CONSEQUENCE = 'major_consequence'
PASS = 'PASS'  # noqa: S105

# Used to filter mt.info fields.
CONSEQ_PREDICTED_PREFIX = 'PREDICTED_'
NON_GENE_PREDICTIONS = {
    'PREDICTED_INTERGENIC',
    'PREDICTED_NONCODING_BREAKPOINT',
    'PREDICTED_NONCODING_SPAN',
}

GENCODE_FILE_HEADER = [
    'chrom',
    'source',
    'feature_type',
    'start',
    'end',
    'score',
    'strand',
    'phase',
    'info',
]


# yoinking some methods out of hail_scripts.computed_fields
# removes dependency on submodule completely
def get_expr_for_contig_number(locus: hl.LocusExpression) -> hl.Int32Expression:
    """Convert contig name to contig number"""
    return hl.bind(
        lambda contig: (
            hl.case().when(contig == 'X', 23).when(contig == 'Y', 24).when(contig[0] == 'M', 25).default(hl.int(contig))
        ),
        locus.contig.replace('^chr', ''),
    )


def get_expr_for_xpos(
    locus: hl.LocusExpression | hl.StructExpression,
) -> hl.Int64Expression:
    """Genomic position represented as a single number = contig_number * 10**9 + position.
    This represents chrom:pos more compactly and allows for easier sorting.
    """
    contig_number = get_expr_for_contig_number(locus)
    return hl.int64(contig_number) * 1_000_000_000 + locus.position


def get_cpx_interval(x):
    """
    an example format of CPX_INTERVALS is "DUP_chr1:1499897-1499974"
    Args:
        x (hl.StringExpression instance):
    Returns:
        Struct of CPX components
    """
    type_chr = x.split('_chr')
    chr_pos = type_chr[1].split(':')
    pos = chr_pos[1].split('-')
    return hl.struct(type=type_chr[0], chrom=chr_pos[0], start=hl.int32(pos[0]), end=hl.int32(pos[1]))


def parse_gtf_from_local(gtf_path: str) -> hl.dict:
    """
    Read over the localised GTF file and read into a dict

    Args:
        gtf_path ():
    Returns:
        the gene lookup dictionary as a Hail DictExpression
    """
    gene_id_mapping = {}
    loguru.logger.info(f'Loading {gtf_path}')
    with gzip.open(gtf_path, 'rt') as gencode_file:
        # iterate over this file and do all the things
        for i, line in enumerate(gencode_file):
            tidy_line = line.rstrip('\r\n')

            if not tidy_line or tidy_line.startswith('#'):
                continue

            fields = tidy_line.split('\t')
            if len(fields) != len(GENCODE_FILE_HEADER):
                raise ValueError(f'Unexpected number of fields on line #{i}: {fields}')
            record = dict(zip(GENCODE_FILE_HEADER, fields, strict=False))
            if record['feature_type'] != 'gene':
                continue
            # parse info field
            info_fields_list = [x.strip().split() for x in record['info'].split(';') if x != '']
            info_fields = {k: v.strip('"') for k, v in info_fields_list}

            # skip an ENSG: ENSG mapping, redundant...
            if info_fields['gene_name'].startswith('ENSG'):
                continue
            gene_id_mapping[info_fields['gene_name']] = info_fields['gene_id'].split('.')[0]
    return hl.literal(gene_id_mapping)


def annotate_cohort_sv(vcf_path: str, out_mt_path: str, gencode_gz: str, checkpoint_path: str):
    """
    Translate an annotated SV VCF into a Seqr-ready format
    Relevant gCNV specific schema
    https://github.com/populationgenomics/seqr-loading-pipelines/blob/master/luigi_pipeline/lib/model/sv_mt_schema.py
    Relevant gCNV loader script
    https://github.com/populationgenomics/seqr-loading-pipelines/blob/master/luigi_pipeline/seqr_sv_loading.py
    Args:
        vcf_path (str): Where is the VCF??
        out_mt_path (str): And where do you need output!?
        gencode_gz (str): The path to a compressed GENCODE GTF file
        checkpoint_path (str): CHECKPOINT!@!!
    """

    hail_batch.init_batch()

    loguru.logger.info(f'Importing SV VCF {vcf_path}')
    mt = hl.import_vcf(
        vcf_path,
        reference_genome=hail_batch.genome_build(),
        skip_invalid_loci=True,
        force_bgz=True,
    )

    # add attributes required for Seqr
    mt = mt.annotate_globals(
        sourceFilePath=vcf_path,
        genomeVersion=hail_batch.genome_build().replace('GRCh', ''),
        hail_version=hl.version(),
        datasetType='SV',
    )

    mt = mt.checkpoint(checkpoint_path)

    if sequencing_type := config.config_retrieve(['workflow', 'sequencing_type']):
        # Map to Seqr-style string
        # https://github.com/broadinstitute/seqr/blob/e0c179c36c0f68c892017de5eab2e4c1b9ffdc92/seqr/models.py#L592-L594
        mt = mt.annotate_globals(
            sampleType={
                'genome': 'WGS',
                'exome': 'WES',
                'single_cell': 'RNA',
            }.get(sequencing_type, ''),
        )

    # reimplementation of
    # github.com/populationgenomics/seqr-loading-pipelines..luigi_pipeline/lib/model/sv_mt_schema.py
    population_prefix = config.config_retrieve(['references', 'gatk_sv', 'external_af_ref_bed_prefix'])
    mt = mt.annotate_rows(
        sc=mt.info.AC[0],
        sf=mt.info.AF[0],
        sn=mt.info.AN,
        end=mt.info.END,
        end_locus=hl.if_else(
            hl.is_defined(mt.info.END2),
            hl.struct(contig=mt.info.CHR2, position=mt.info.END2),
            hl.struct(contig=mt.locus.contig, position=mt.info.END),
        ),
        sv_callset_Het=mt.info.N_HET,
        sv_callset_Hom=mt.info.N_HOMALT,
        gnomad_svs_ID=mt.info[f'{population_prefix}_SVID'],
        gnomad_svs_AF=mt.info[f'{population_prefix}_AF'],
        gnomad_svs_AC=hl.missing('float64'),
        gnomad_svs_AN=hl.missing('float64'),
        StrVCTVRE_score=hl.parse_float(mt.info.StrVCTVRE),
        filters=hl.or_missing(
            (mt.filters.filter(lambda x: (x != PASS) & (x != BOTHSIDES_SUPPORT))).size() > 0,
            mt.filters,
        ),
        bothsides_support=mt.filters.any(lambda x: x == BOTHSIDES_SUPPORT),
        algorithms=mt.info.ALGORITHMS,
        cpx_intervals=hl.or_missing(
            hl.is_defined(mt.info.CPX_INTERVALS),
            mt.info.CPX_INTERVALS.map(lambda x: get_cpx_interval(x)),
        ),
        sv_types=mt.alleles[1].replace('[<>]', '').split(':', 2),
    )

    # get the Gene-Symbol mapping dict
    gene_id_mapping = parse_gtf_from_local(gencode_gz)[0]

    # OK, NOW IT'S BUSINESS TIME
    conseq_predicted_gene_cols = [
        gene_col
        for gene_col in mt.info
        if (gene_col.startswith(CONSEQ_PREDICTED_PREFIX) and gene_col not in NON_GENE_PREDICTIONS)
    ]

    # register a chain file
    liftover_path = config.config_retrieve(['references', 'liftover_38_to_37'])
    rg37 = hl.get_reference('GRCh37')
    rg38 = hl.get_reference('GRCh38')
    rg38.add_liftover(str(liftover_path), rg37)

    # annotate with mapped genes
    # Note I'm adding a Flake8 noqa for B023 (loop variable gene_col unbound)
    # I've experimented in a notebook and this seems to perform as expected
    # The homologous small variant seqr_loader method performs a similar function
    # but in a slightly more complicated way (mediated by a method in the S-L-P
    # codebase, so as not to trigger Flake8 evaluation)
    # pos/contig/xpos sourced from
    # seqr-loading-pipelines...luigi_pipeline/lib/model/seqr_mt_schema.py#L12
    mt = mt.annotate_rows(
        sortedTranscriptConsequences=hl.filter(
            hl.is_defined,
            [
                mt.info[gene_col].map(
                    lambda gene: hl.struct(
                        **{
                            GENE_SYMBOL: gene,
                            GENE_ID: gene_id_mapping.get(gene, hl.missing(hl.tstr)),
                            MAJOR_CONSEQUENCE: gene_col.replace(CONSEQ_PREDICTED_PREFIX, '', 1),  # noqa: B023
                        },
                    ),
                )
                for gene_col in conseq_predicted_gene_cols
            ],
        ).flatmap(lambda x: x),
        contig=mt.locus.contig.replace('^chr', ''),
        start=mt.locus.position,
        pos=mt.locus.position,
        xpos=get_expr_for_xpos(mt.locus),
        xstart=get_expr_for_xpos(mt.locus),
        xstop=get_expr_for_xpos(mt.end_locus),
        rg37_locus=hl.liftover(mt.locus, 'GRCh37'),
        rg37_locus_end=hl.or_missing(
            mt.end_locus.position <= hl.literal(hl.get_reference('GRCh38').lengths)[mt.end_locus.contig],
            hl.liftover(
                hl.locus(
                    mt.end_locus.contig,
                    mt.end_locus.position,
                    reference_genome='GRCh38',
                ),
                'GRCh37',
            ),
        ),
        svType=mt.sv_types[0],
        sv_type_detail=hl.if_else(
            mt.sv_types[0] == 'CPX',
            mt.info.CPX_TYPE,
            hl.or_missing((mt.sv_types[0] == 'INS') & (hl.len(mt.sv_types) > 1), mt.sv_types[1]),
        ),
        variantId=mt.rsid,
        docId=mt.rsid[0:512],
    )

    # and some more annotation stuff
    mt = mt.annotate_rows(
        transcriptConsequenceTerms=hl.set(
            mt.sortedTranscriptConsequences.map(lambda x: x[MAJOR_CONSEQUENCE]).extend([mt.sv_types[0]]),
        ),
        geneIds=hl.set(
            mt.sortedTranscriptConsequences.filter(lambda x: x[MAJOR_CONSEQUENCE] != 'NEAREST_TSS').map(
                lambda x: x[GENE_ID],
            ),
        ),
        rsid=hl.missing('tstr'),
    )

    # write this output
    mt.write(out_mt_path, overwrite=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Annotate cohort SV VCF')
    parser.add_argument('--vcf', required=True, help='Path to the input VCF file')
    parser.add_argument('--output', required=True, help='Path to the output MT file')
    parser.add_argument('--gencode_gtf', required=True, help='Path to the GENCODE GTF file (gzipped)')
    parser.add_argument('--checkpoint', required=True, help='Checkpoint for MT')

    args = parser.parse_args()

    annotate_cohort_sv(args.vcf, args.output, args.gencode_gtf, args.checkpoint)
