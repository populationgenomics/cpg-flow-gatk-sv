from argparse import ArgumentParser

from cpg_utils import config, metamist_registration, to_path

def read_file_lines(filepath: str) -> list[str]:
    """Read contents of a single file."""
    contents: list[str] = []
    with to_path(filepath).open() as handle:
        for line in handle:
            contents.append(line.strip())
    return contents

def main(
        dataset: str,
        output: str,
        stage: str,
        analysis_type: str,
        exclusion_file: str,
        sg_file: str,
):
    exclusions = read_file_lines(exclusion_file)
    sgs = read_file_lines(sg_file)

    included = [sg for sg in sgs if sg not in exclusions]
    relevant_exclusions = [sg for sg in exclusions if sg in sgs]

    metamist_registration.create_new(
        project=config.dataset_for_access_level(dataset),
        output=output,
        analysis_type=analysis_type,
        meta={'stage': stage, 'sequencing_type': 'genome', 'exclusions': relevant_exclusions},
        sgs=included,
    )



if __name__ == '__main__':
    parser = ArgumentParser()
    parser.add_argument('--dataset', help='Project to register the data into.', required=True)
    parser.add_argument('--output', help='Path to Stage output.', required=True)
    parser.add_argument('--stage', help='Stage the analysis originates from.', required=True)
    parser.add_argument('--atype', help='Analysis type to create.', required=True)
    parser.add_argument('--sgs',help='Path to a file containing SG IDs.', required=True)
    parser.add_argument('--exclusions', help='Path to excluded-SGs file.', required=True)
    args = parser.parse_args()

