"""
Hail Query functions for seqr loader; SV edition.
"""

import gzip
import loguru

import hail as hl

from cpg_workflows.utils import generator_chunks, read_hail


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


def annotate_dataset_sv(mt_path: str, out_mt_path: str):
    """
    load the stuff specific to samples in this dataset
    do this after subsetting to specific samples

    Removing the current logic around comparing genotypes to a previous
    callset - doesn't fit with the current implementation

    Args:
        mt_path (str): path to the annotated MatrixTable
        out_mt_path (str): and where do you want it to end up?
    """

    loguru.logger.info('Annotating genotypes')

    mt = read_hail(mt_path)
    is_called = hl.is_defined(mt.GT)
    num_alt = hl.if_else(is_called, mt.GT.n_alt_alleles(), -1)
    mt = mt.annotate_rows(
        genotypes=hl.agg.collect(
            hl.struct(
                sample_id=mt.s,
                gq=mt.GQ,
                cn=mt.RD_CN,
                num_alt=num_alt,
                # prev_num_alt=hl.or_missing(discordant_genotype, prev_num_alt),
                # prev_call=hl.or_missing(
                #     is_called, was_previously_called & concordant_genotype
                # ),
                # new_call=hl.or_missing(
                #     is_called, ~was_previously_called | novel_genotype
                # ),
            ),
        ),
    )

    def _genotype_filter_samples(fn):
        # Filter on the genotypes.
        return hl.set(mt.genotypes.filter(fn).map(lambda g: g.sample_id))

    # top level - decorator
    def _capture_i_decorator(func):
        # call the returned_function(i) which locks in the value of i
        def _inner_filter(i):
            # the _genotype_filter_samples will call this _func with g
            def _func(g):
                return func(i, g)

            return _func

        return _inner_filter

    @_capture_i_decorator
    def _filter_num_alt(i, g):
        return i == g.num_alt

    @_capture_i_decorator
    def _filter_samples_gq(i, g):
        return (g.gq >= i) & (g.gq < i + 10)

    @_capture_i_decorator
    def _filter_sample_cn(i, g):
        return g.cn == i

    # github.com/populationgenomics/seqr-loading-pipelines/blob/master/luigi_pipeline/lib/model/sv_mt_schema.py#L221
    # github.com/populationgenomics/seqr-loading-pipelines/blob/master/luigi_pipeline/lib/model/seqr_mt_schema.py#L251
    mt = mt.annotate_rows(
        samples_no_call=_genotype_filter_samples(lambda g: g.num_alt == -1),
        samples_num_alt=hl.struct(**{'%i' % i: _genotype_filter_samples(_filter_num_alt(i)) for i in range(1, 3, 1)}),
        samples_gq_sv=hl.struct(
            **{('%i_to_%i' % (i, i + 10)): _genotype_filter_samples(_filter_samples_gq(i)) for i in range(0, 90, 10)},
        ),
    )

    loguru.logger.info('Genotype fields annotated')
    mt.describe()
    mt.write(out_mt_path, overwrite=True)
    loguru.logger.info(f'Written  SV MT to {out_mt_path}')


if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser(description='Annotate SV dataset MatrixTable')
    parser.add_argument('--mt', required=True, help='Path to the input MatrixTable')
    parser.add_argument('--out_mt', required=True, help='Path to write the annotated MatrixTable')
    args = parser.parse_args()

    annotate_dataset_sv(args.mt, args.out_mt)
