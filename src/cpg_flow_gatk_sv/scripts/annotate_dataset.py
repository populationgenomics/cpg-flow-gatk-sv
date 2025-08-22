"""
Hail Query functions for seqr loader; SV edition.
"""

import argparse

import loguru

import hail as hl

from cpg_utils import hail_batch


def subset_mt_to_samples(mt_path: str, sample_id_file: str, dataset_mt: str, exclusion_file: str | None = None):
    """
    Subset the MatrixTable to the provided list of samples and to variants present
    in those samples

    Args:
        mt_path (str): cohort-level matrix table from VCF.
        sample_id_file (list[str]): samples to take from the matrix table, one per line
        dataset_mt (str): path to write the result.
        exclusion_file (str, optional): path to a file containing samples to remove from the
                                        subset prior to extracting
    """

    loguru.logger.info(f'Subsetting cohort MatrixTable to samples, reading from {mt_path}')
    mt = hl.read_matrix_table(mt_path)

    # read the sample IDs from the file
    with hl.hadoop_open(sample_id_file) as f:
        sample_ids = set(f.read().splitlines())

    # if an exclusion file was passed, remove the samples from the subset
    # this executes in a query command, by execution time the file should exist
    if exclusion_file:
        with hl.hadoop_open(exclusion_file) as f:
            exclusion_ids = set(f.read().splitlines())
        loguru.logger.info(f'Excluding {len(exclusion_ids)} samples from the subset')
        sample_ids -= exclusion_ids

    mt_sample_ids = set(mt.s.collect())

    if sample_ids_not_in_mt := sample_ids - mt_sample_ids:
        raise ValueError(
            f'Found {len(sample_ids_not_in_mt)}/{len(sample_ids)} IDs in the requested subset not in the callset.\n'
            f"IDs that aren't in the callset: {sample_ids_not_in_mt}\n"
            f'All callset sample IDs: {mt_sample_ids}',
        )

    loguru.logger.info(f'Found {len(mt_sample_ids)} samples in mt, subsetting to {len(sample_ids)} samples.')

    mt = mt.filter_cols(hl.literal(sample_ids).contains(mt.s))
    mt = mt.filter_rows(hl.agg.any(mt.GT.is_non_ref()))
    mt.write(dataset_mt, overwrite=True)
    loguru.logger.info(f'Finished subsetting to {len(sample_ids)} samples, written to {dataset_mt}.')


def annotate_dataset(mt_path: str, out_mt_path: str):
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

    mt = hl.read_matrix_table(mt_path)
    is_called = hl.is_defined(mt.GT)
    num_alt = hl.if_else(is_called, mt.GT.n_alt_alleles(), -1)
    mt = mt.annotate_rows(
        genotypes=hl.agg.collect(
            hl.struct(
                sample_id=mt.s,
                gq=mt.GQ,
                cn=mt.RD_CN,
                num_alt=num_alt,
            ),
        ),
    )

    def _genotype_filter_samples(fn) -> hl.set:
        # Filter on the genotypes.
        return hl.set(mt.genotypes.filter(fn).map(lambda g: g.sample_id))

    # top level - decorator
    def _capture_i_decorator(func):  # noqa: ANN202
        # call the returned_function(i) which locks in the value of i
        def _inner_filter(i):  # noqa: ANN202
            # the _genotype_filter_samples will call this _func with g
            def _func(g):  # noqa: ANN202
                return func(i, g)

            return _func

        return _inner_filter

    @_capture_i_decorator
    def _filter_num_alt(i, g) -> hl.bool:
        return i == g.num_alt

    @_capture_i_decorator
    def _filter_samples_gq(i, g) -> hl.bool:
        return (g.gq >= i) & (g.gq < i + 10)

    @_capture_i_decorator
    def _filter_sample_cn(i, g) -> hl.bool:
        return g.cn == i

    # github.com/populationgenomics/seqr-loading-pipelines/blob/master/luigi_pipeline/lib/model/sv_mt_schema.py#L221
    # github.com/populationgenomics/seqr-loading-pipelines/blob/master/luigi_pipeline/lib/model/seqr_mt_schema.py#L251
    mt = mt.annotate_rows(
        samples_no_call=_genotype_filter_samples(lambda g: g.num_alt == -1),
        samples_num_alt=hl.struct(**{f'{i}': _genotype_filter_samples(_filter_num_alt(i)) for i in range(1, 3, 1)}),
        samples_gq_sv=hl.struct(
            **{f'{i}_to_{i + 10}': _genotype_filter_samples(_filter_samples_gq(i)) for i in range(0, 90, 10)},
        ),
    )

    loguru.logger.info('Genotype fields annotated')
    mt.describe()
    mt.write(out_mt_path, overwrite=True)
    loguru.logger.info(f'Written  SV MT to {out_mt_path}')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Annotate SV dataset MatrixTable')
    parser.add_argument('--mt', required=True, help='Path to the input MatrixTable')
    parser.add_argument('--dataset_mt', required=True, help='Path to write the single-dataset MT')
    parser.add_argument('--output', required=True, help='Path to write the final MatrixTable')
    parser.add_argument(
        '--sample_id_file',
        required=True,
        help='Path to a file containing sample IDs to include in the subset, one per line',
    )
    parser.add_argument(
        '--exclusion_file',
        required=False,
        help='Path to a file containing sample IDs to exclude from the subset, one per line',
    )
    args = parser.parse_args()

    hail_batch.init_batch()

    subset_mt_to_samples(
        mt_path=args.mt,
        sample_id_file=args.sample_id_file,
        dataset_mt=args.dataset_mt,
        exclusion_file=args.exclusion_file,
    )

    annotate_dataset(
        args.dataset_mt,
        args.output,
    )
