from typing import TYPE_CHECKING

from cpg_flow import targets
from cpg_utils import Path, config, hail_batch

if TYPE_CHECKING:
    from hailtop.batch.job import BashJob


def create_strip_qname_suffix_job(
    sg: targets.SequencingGroup,
    expected_outputs: dict[str, Path],
) -> 'BashJob':
    """
    Strip /1 and /2 QNAME suffixes from a CRAM, preserving coordinate sort order.

    Some older sequencing runs leave /1 /2 suffixes on read QNAMEs. These break
    mate pairing in any downstream step that runs samtools fastq without collation
    (e.g. RealignSoftClippedReads in the GATK-SV WDL).
    """
    batch = hail_batch.get_batch()
    j = batch.new_bash_job(f'StripQnameSuffixes_{sg.id}')
    j.image(config.config_retrieve(['images', 'samtools_cloud_docker']))
    j.cpu(4)
    j.memory('standard')
    j.storage('80Gi')

    reference = hail_batch.fasta_res_group(batch)

    input_cram = batch.read_input_group(
        cram=str(sg.cram),
        **{'cram.crai': str(sg.cram) + '.crai'},
    )

    output_cram = str(expected_outputs['cram'])
    output_crai = str(expected_outputs['crai'])

    awk_strip = r"""awk 'BEGIN{FS=OFS="\t"} !/^@/{sub(/\/[12]$/,"",$1)} {print}'"""

    cmd = f"""\
samtools view -h -T {reference.base} -@ 3 {input_cram.cram} | \
{awk_strip} | \
samtools view -C -T {reference.base} -@ 3 -o /tmp/cleaned.cram -

samtools index /tmp/cleaned.cram /tmp/cleaned.cram.crai

gcloud storage cp /tmp/cleaned.cram {output_cram}
gcloud storage cp /tmp/cleaned.cram.crai {output_crai}
"""

    j.command(hail_batch.command(cmd, setup_gcp=True))
    return j