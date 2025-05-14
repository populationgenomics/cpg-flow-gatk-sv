#!/usr/bin/env python3

from argparse import ArgumentParser

from cpg_flow.workflow import run_workflow

from cpg_flow_gatk_sv.first_stages import GatherSampleEvidence


def cli_main():
    """
    CLI entrypoint - starts up the workflow
    """
    parser = ArgumentParser()
    parser.add_argument('--dry_run', action='store_true', help='Dry run')
    args = parser.parse_args()

    stages = [GatherSampleEvidence]

    run_workflow(stages=stages, dry_run=args.dry_run)


if __name__ == '__main__':
    cli_main()
