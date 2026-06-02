# tileylab/pgp: Usage

> _Documentation of pipeline parameters is generated automatically from the pipeline schema and can no longer be found in markdown files._

## Introduction

This page documents the samplesheet format, the typical command line, and the available profiles. Parameter reference is generated automatically from the schema — run `nextflow run tileylab/pgp --help` for the full list.

## Samplesheet input

You will need to create a samplesheet describing your samples before running the pipeline. It must be a comma-separated file with a header row and five columns:

```bash
--input '[path to samplesheet file]'
```

| Column    | Description                                                                                       |
| --------- | ------------------------------------------------------------------------------------------------- |
| `sample`  | Unique sample identifier (no spaces). Becomes `meta.id` in the workflow.                          |
| `fastq_1` | Path to gzipped FastQ for read 1 (`.fastq.gz` or `.fq.gz`).                                       |
| `fastq_2` | Path to gzipped FastQ for read 2. Leave empty for single-end data.                                |
| `species` | Species, population, or other grouping identifier (no spaces). Required.                          |
| `ploidy`  | Sample ploidy as a positive integer. Required.                                                    |

`species` and `ploidy` are both required by [`assets/schema_input.json`](../assets/schema_input.json); runs will fail validation if either is missing.

### Multiple lanes for the same sample

The same `sample` identifier may appear on multiple rows when a library was sequenced across lanes. The pipeline groups by `sample` and concatenates raw reads before downstream analysis. All rows for a given sample must agree on endedness (all paired or all single).

```csv title="samplesheet.csv"
sample,fastq_1,fastq_2,species,ploidy
CONTROL_REP1,AEG588A1_S1_L002_R1_001.fastq.gz,AEG588A1_S1_L002_R2_001.fastq.gz,Vdarrowii,2
CONTROL_REP1,AEG588A1_S1_L003_R1_001.fastq.gz,AEG588A1_S1_L003_R2_001.fastq.gz,Vdarrowii,2
CONTROL_REP1,AEG588A1_S1_L004_R1_001.fastq.gz,AEG588A1_S1_L004_R2_001.fastq.gz,Vdarrowii,2
```

### Mixed single- and paired-end example

```csv title="samplesheet.csv"
sample,fastq_1,fastq_2,species,ploidy
CONTROL_REP1,AEG588A1_S1_L002_R1_001.fastq.gz,AEG588A1_S1_L002_R2_001.fastq.gz,Vdarrowii,2
CONTROL_REP2,AEG588A2_S2_L002_R1_001.fastq.gz,AEG588A2_S2_L002_R2_001.fastq.gz,Vdarrowii,2
TREATMENT_REP1,AEG588A4_S4_L003_R1_001.fastq.gz,,Vcorymbosum,4
TREATMENT_REP2,AEG588A5_S5_L003_R1_001.fastq.gz,,Vcorymbosum,4
```

An [example samplesheet](../assets/samplesheet.csv) ships with the pipeline.

## Running the pipeline

A typical command:

```bash
nextflow run tileylab/pgp \
    --input ./samplesheet.csv \
    --outdir ./results \
    --reference ./ref/ref.fa \
    --intervals ./ref/intervals.list \
    -profile docker
```

`--reference` and `--intervals` are required in addition to `--input` and `--outdir`. The PhiX contamination reference used by `BBMAP_BBDUK` ships with the pipeline at `assets/phix174_ill.ref.fa` and is not configurable on the command line.

This will launch the pipeline with the `docker` configuration profile. See below for more profiles.

Note that the pipeline will create the following files in your working directory:

```bash
work                # Directory containing the nextflow working files
<OUTDIR>            # Finished results in specified location (defined with --outdir)
.nextflow_log       # Log file from Nextflow
# Other nextflow hidden files, eg. history of pipeline runs and old logs.
```

If you wish to repeatedly use the same parameters for multiple runs, rather than specifying each flag in the command, you can specify these in a params file. Pipeline settings can be provided in a `yaml` or `json` file via `-params-file <file>`.

> [!WARNING]
> Do not use `-c <file>` to specify parameters as this will result in errors. Custom config files specified with `-c` must only be used for [tuning process resource specifications](https://www.nextflow.io/docs/latest/config.html#scope-process), other infrastructural tweaks (such as output directories), or module arguments (args).

The above pipeline run specified with a params file in yaml format:

```bash
nextflow run tileylab/pgp -profile docker -params-file params.yaml
```

with:

```yaml title="params.yaml"
input: './samplesheet.csv'
outdir: './results/'
reference: './ref/ref.fa'
intervals: './ref/intervals.list'
```

### Updating the pipeline

When you run the above command, Nextflow automatically pulls the pipeline code from GitHub and stores it as a cached version. When running the pipeline after this, it will always use the cached version if available - even if the pipeline has been updated since. To make sure that you're running the latest version of the pipeline, make sure that you regularly update the cached version of the pipeline:

```bash
nextflow pull tileylab/pgp
```

### Reproducibility

It is a good idea to specify the pipeline version when running the pipeline on your data. This ensures that a specific version of the pipeline code and software are used when you run your pipeline. If you keep using the same tag, you'll be running the same version of the pipeline, even if there have been changes to the code since.

First, go to the [tileylab/pgp releases page](https://github.com/tileylab/pgp/releases) and find the latest pipeline version - numeric only (eg. `1.3.1`). Then specify this when running the pipeline with `-r` (one hyphen) - eg. `-r 1.3.1`. Of course, you can switch to another version by changing the number after the `-r` flag.

This version number will be logged in reports when you run the pipeline, so that you'll know what you used when you look back in the future. For example, at the bottom of the MultiQC reports.

To further assist in reproducibility, you can use share and reuse [parameter files](#running-the-pipeline) to repeat pipeline runs with the same settings without having to write out a command with every single parameter.

## Core Nextflow arguments

> [!NOTE]
> These options are part of Nextflow and use a _single_ hyphen (pipeline parameters use a double-hyphen).

### `-profile`

Use this parameter to choose one or more configuration profiles. Profiles can be combined (e.g. `-profile slurm,singularity`); later profiles override earlier ones.

> [!IMPORTANT]
> We highly recommend Docker or Singularity containers for full pipeline reproducibility. Conda is supported as a fallback.

Profiles bundled with the pipeline:

**Containers / package managers**

- `docker` — [Docker](https://docker.com/)
- `singularity` — [Singularity](https://sylabs.io/docs/)
- `apptainer` — [Apptainer](https://apptainer.org/)
- `podman` — [Podman](https://podman.io/)
- `shifter` — [Shifter](https://nersc.gitlab.io/development/shifter/how-to-use/)
- `charliecloud` — [Charliecloud](https://charliecloud.io/)
- `wave` — [Wave](https://seqera.io/wave/) dynamic container builds (combine with another container profile)
- `conda` / `mamba` — [Conda](https://conda.io/docs/) / [Mamba](https://mamba.readthedocs.io/) (use as a last resort)

**Executors**

- `slurm` — submit to a Slurm cluster. Pin the partition with `--slurm_queue <partition>`; otherwise the cluster default queue is used. `executor.queueSize` is capped at 100 and submissions are rate-limited to 10/min.

**Architecture / accelerators**

- `arm64` — build/pull arm64 containers via Wave
- `emulate_amd64` — run amd64 containers under Docker's arm64 emulation
- `gpu` — pass GPU flags to docker/singularity/apptainer

**Testing**

- `test` — runs the synthetic dataset under `tests/data/`

If `-profile` is not specified, the pipeline will run locally and expect all software to be installed and available on the `PATH`. This is _not_ recommended.

### `-resume`

Specify this when restarting a pipeline. Nextflow will use cached results from any pipeline steps where the inputs are the same, continuing from where it got to previously. For input to be considered the same, not only the names must be identical but the files' contents as well. For more info about this parameter, see [this blog post](https://www.nextflow.io/blog/2019/demystifying-nextflow-resume.html).

You can also supply a run name to resume a specific run: `-resume [run-name]`. Use the `nextflow log` command to show previous run names.

### `-c`

Specify the path to a specific config file (this is a core Nextflow command). Useful for institutional resource tweaks (e.g. a custom `slurm.config` with per-process queues, account strings, or scratch directories).

## Custom configuration

### Resource requests

Whilst the default requirements set within the pipeline will hopefully work for most people and with most input data, you may find that you want to customise the compute resources that the pipeline requests. Each step in the pipeline has a default set of requirements for number of CPUs, memory, and time. Per-process overrides live in `conf/modules.config`; default per-label resources live in `conf/base.config`. Failed jobs are automatically resubmitted with higher resources (2x then 3x) on retryable exit codes.

### Custom containers

To use a different container or conda environment than the one specified in a module, override it in your own `-c` config using a process selector. See the [Nextflow process selector docs](https://www.nextflow.io/docs/latest/config.html#config-process-selectors).

### Custom tool arguments

Each module exposes `ext.args` (and `ext.args2`/`ext.args3` where applicable). Override in a custom `-c` config:

```groovy
process {
    withName: 'BWA_MEM' {
        ext.args = '-M -t 8'
    }
}
```

## Running in the background

Nextflow handles job submissions and supervises the running jobs. The Nextflow process must run until the pipeline is finished.

The Nextflow `-bg` flag launches Nextflow in the background, detached from your terminal so that the workflow does not stop if you log out of your session. The logs are saved to a file.

Alternatively, you can use `screen` / `tmux` or similar tool to create a detached session which you can log back into at a later time. Some HPC setups also allow you to run nextflow within a cluster job submitted your job scheduler (from where it submits more jobs).

## Nextflow memory requirements

In some cases, the Nextflow Java virtual machines can start to request a large amount of memory. We recommend adding the following line to your environment to limit this (typically in `~/.bashrc` or `~./bash_profile`):

```bash
NXF_OPTS='-Xms1g -Xmx4g'
```
