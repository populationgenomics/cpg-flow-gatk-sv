# GATK-SV, in CPG-Flow

Current Version: 0.1.24

This repository contains the GATK-SV workflow, which is a wrapper around the [GATK-SV Cromwell workflow](https://github.com/broadinstitute/gatk-sv). It has been migrated from [Production-Pipelines](https://github.com/populationgenomics/production-pipelines/tree/main/cpg_workflows/stages/gatk_sv). This has been generated as part of the migration from Production-Pipelines's CPG_workflows framework, to the separate CPG-Flow.

## Workflows

This repository contains the following workflows:

- `singlesample_workflow` - previously called `gatk_sv_singlesample`. This is the initial per-sample variant calling, generation of QC metrics following variant calling, through to batching of samples based on those metrics. The exact stages are:
    1. [GatherSampleEvidence](https://github.com/broadinstitute/gatk-sv/blob/main/wdl/GatherSampleEvidence.wdl)
    2. [EvidenceQC](https://github.com/broadinstitute/gatk-sv/blob/main/wdl/EvidenceQC.wdl)
    3. `CreateSampleBatches` - a somewhat bespoke sample batching script

- `multisample_workflow` - previously called `gatk_sv_multisample`. This is the Cohort based batch processing of clustered samples, through to a final joint-called dataset.

## Structure

Following a best-practices CPG-Flow structure, the repository is structured as follows:

```commandline
├── Dockerfile
├── LICENSE
├── README.md
├── pull_request_template.md
├── pyproject.toml
└── src
    ├── cpg_flow_gatk_sv
    │   ├── __init__.py
    │   ├── config_template.toml
    │   ├── singlesample_workflow.py
    │   ├── multisample_workflow.py
    │   ├── jobs
    │   │   ├── CreateSampleBatches.py
    │   │   ├── EvidenceQC.py
    │   │   ├── GatherSampleEvidence.py
    │   │   └── ...
    │   ├── scripts
    │   │   ├── __init__.py
    │   │   └── sample_batching.py
    │   │   └── ...
    │   └── utils.py
```

This structure contains the following important files:

- `pyproject.toml` - the build instructions and linter settings for the repository
- `singlesample_workflow` - the Stages and workflow trigger for the first workflow
- `multisample_workflow` - the Stages and workflow trigger for the second workflow
- `jobs` - a directory containing the logic for each of the Stages in the workflow. Each file contains a single Stage, and is
  named after the Stage it contains. This implements all the logic for a stage, including assembling the input dictionary to be passed to Cromwell.
- `scripts` - a directory containing any scripts which are not part of the GATK-SV Cromwell workflow, but are used in the workflow. This includes any helper scripts, or CPG-custom logic.

## Example invocation

This is designed to be run using analysis-runner, using a fully containerised installation of this codebase.

```bash
analysis-runner \
    --skip-repo-checkout \
    --image australia-southeast1-docker.pkg.dev/cpg-common/images/cpg-flow-gatk-sv:0.1.24 \
    --dataset DATASET \
    --description 'GATK-SV, CPG-flow' \
    -o gatk-sv_cpg-flow \
    --access-level full \
    --config CONFIG \
    singlesample_workflow
```

## Development

Semantic Versioning should be implemented with `bump-my-version`

```commandline
bump-my-version bump patch/minor/major
```

This is using the configuration block inside `pyproject.toml`
