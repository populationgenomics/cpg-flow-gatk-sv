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

    subprocess.run(
        ['gcloud', 'auth', 'configure-docker', 'australia-southeast1-docker.pkg.dev'],
        check=True,
    )

    response = urlopen(args.docker_json)
    dockers_json = json.loads(response.read())
    dockers_json.pop('name')

    config_section = {}
    for i, key in enumerate(dockers_json):
        image_name = dockers_json[key].split('/')[-1]
        cpg_ar_path = 'australia-southeast1-docker.pkg.dev/cpg-common/images/sv/' + image_name
        print(f'#{i}: copying {key}: {dockers_json[key]} to {cpg_ar_path}')
        src_path = 'docker://' + dockers_json[key]
        dst_path = 'docker://' + cpg_ar_path
        cmd = f'skopeo copy {src_path} {dst_path}'
        subprocess.run(cmd, shell=True, check=True)
        config_section[key] = cpg_ar_path

    print()
    print('TOML [images] config section:')
    print()
    print(toml.dumps(config_section))


if __name__ == '__main__':
    main()
