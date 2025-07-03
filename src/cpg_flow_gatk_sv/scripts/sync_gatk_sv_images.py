"""
A simple script to copy docker images specified in a json file to CPG artefact registry.

Typical usage has been pointing this script to the raw content of a publicly available JSON file
e.g. https://raw.githubusercontent.com/populationgenomics/gatk-sv/main/inputs/values/dockers.json
"""

import argparse
import json
import subprocess
from urllib.request import urlopen

import toml


def main():
    """Copies each docker image defined in the JSON argument."""
    parser = argparse.ArgumentParser(description='Sync GATK SV docker images to CPG artefact registry.')
    parser.add_argument(
        'docker_json',
        help='URL to the JSON file containing docker images to copy.',
    )
    args = parser.parse_args()

    subprocess.run(  # noqa: S603
        ['gcloud', 'auth', 'configure-docker', 'australia-southeast1-docker.pkg.dev'],  # noqa: S607
        check=True,
    )

    response = urlopen(args.docker_json)  # noqa: S310
    dockers_json = json.loads(response.read())

    config_section = {}
    for key, value in dockers_json.items():
        image_name = value.split('/')[-1]
        cpg_ar_path = 'australia-southeast1-docker.pkg.dev/cpg-common/images/sv/' + image_name
        print(f'Copying {key}: {value} to {cpg_ar_path}')
        subprocess.run(f'skopeo copy docker://{value} docker://{cpg_ar_path}', shell=True, check=True)  # noqa: S602
        config_section[key] = cpg_ar_path

    print()
    print('TOML [images] config section:')
    print()
    print(toml.dumps(config_section))


if __name__ == '__main__':
    main()
