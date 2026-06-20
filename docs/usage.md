# tileylab/pgp: Usage

> _Documentation of pipeline parameters is generated automatically from the pipeline schema and can no longer be found in markdown files._

## Introduction

This page documents the samplesheet format, the typical command line, and the available profiles. Parameter reference is generated automatically from the schema — run `nextflow run tileylab/pgp --help` for the full list.

## Samplesheet input

You will need to create a samplesheet describing your samples before running the pipeline. It must be a comma-separated file with a header row and six columns:

```bash
--input '[path to samplesheet file]'
```

| Column     | Description                                                                                               |
| ---------- | --------------------------------------------------------------------------------------------------------- |
| `sample`   | Unique sample identifier (no spaces). Becomes `meta.id` in the workflow.                                  |
| `fastq_1`  | Path to gzipped FastQ for read 1 (`.fastq.gz` or `.fq.gz`).                                               |
| `fastq_2`  | Path to gzipped FastQ for read 2. Leave empty for single-end data.                                        |
| `species`  | Species, population, or other grouping identifier (no spaces). Required.                                  |
| `ploidy`   | Sample ploidy as a positive integer. Required.                                                            |
| `outgroup` | `0` (ingroup) or `1` (outgroup). Required. Outgroups are removed from the ingroup-only VCFtools products. |

`species`, `ploidy`, and `outgroup` are all required by [`assets/schema_input.json`](../assets/schema_input.json); runs will fail validation if any is missing.

### Multiple lanes for the same sample

The same `sample` identifier may appear on multiple rows when a library was sequenced across lanes. The pipeline groups by `sample` and concatenates raw reads before downstream analysis. All rows for a given sample must agree on endedness (all paired or all single).

```csv title="samplesheet.csv"
sample,fastq_1,fastq_2,species,ploidy,outgroup
CONTROL_REP1,AEG588A1_S1_L002_R1_001.fastq.gz,AEG588A1_S1_L002_R2_001.fastq.gz,Vdarrowii,2,0
CONTROL_REP1,AEG588A1_S1_L003_R1_001.fastq.gz,AEG588A1_S1_L003_R2_001.fastq.gz,Vdarrowii,2,0
CONTROL_REP1,AEG588A1_S1_L004_R1_001.fastq.gz,AEG588A1_S1_L004_R2_001.fastq.gz,Vdarrowii,2,0
```

### Mixed single- and paired-end example

```csv title="samplesheet.csv"
sample,fastq_1,fastq_2,species,ploidy,outgroup
CONTROL_REP1,AEG588A1_S1_L002_R1_001.fastq.gz,AEG588A1_S1_L002_R2_001.fastq.gz,Vdarrowii,2,0
CONTROL_REP2,AEG588A2_S2_L002_R1_001.fastq.gz,AEG588A2_S2_L002_R2_001.fastq.gz,Vdarrowii,2,0
TREATMENT_REP1,AEG588A4_S4_L003_R1_001.fastq.gz,,Vcorymbosum,4,1
TREATMENT_REP2,AEG588A5_S5_L003_R1_001.fastq.gz,,Vcorymbosum,4,1
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
input: "./samplesheet.csv"
outdir: "./results/"
reference: "./ref/ref.fa"
intervals: "./ref/intervals.list"
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

### Skip Duplicate Marking For Targeted Libraries

For assays where duplicate marking is not biologically appropriate (for example restriction-enzyme or probe-targeted sequencing), you can bypass Picard duplicate marking:

```bash
nextflow run tileylab/pgp \
  --input ./samplesheet.csv \
  --outdir ./results \
  --reference ./reference.fa \
  --intervals ./intervals.list \
  --skip_markduplicates \
  -profile docker
```

When `--skip_markduplicates` is enabled, the workflow sends `SAMTOOLS_VIEW`-filtered BAMs directly to HaplotypeCaller after indexing, instead of running `PICARD_MARKDUPLICATES`.

## Hard filtering strategy

Variant filtering runs in two stages after joint genotyping. Both stages are driven entirely by `ext.args` in [`conf/modules.config`](../conf/modules.config), so any threshold or flag can be changed in one place (see [Tuning filters](#tuning-filters-via-extargs) below).

### Stage 1 — GATK best-practice hard filters

The multi-sample VCF from `GATK4_GENOTYPEGVCFS` is filtered following GATK best practices:

- `SELECT_ALL` keeps SNPs only (`--select-type-to-include SNP --exclude-filtered`) → `*.snps.vcf.gz`
- `GATK4_VARIANTFILTRATION` annotates the `FILTER` column for sites failing any threshold below (nothing is dropped yet) → `*.snps.annotations.vcf.gz`
- `SELECT_PASS` keeps only `PASS` sites (`--select-type-to-include SNP --exclude-filtered`) → `*.snps.gatkfilters.vcf.gz` — **this is the input to the VCFtools stage**
- `SELECT_BIPASS` additionally restricts to biallelic SNPs → `*.biallelic.gatkfilters.vcf.gz`

| Filter expression       | Flag name             | Annotation               |
| ----------------------- | --------------------- | ------------------------ |
| `FS > 60.0`             | `FS_gt60`             | Fisher strand bias       |
| `MQ < 40.0`             | `MQ_lt40`             | RMS mapping quality      |
| `ReadPosRankSum < -8.0` | `ReadPosRankSum_ltm8` | read-position rank-sum   |
| `MQRankSum < -12.5`     | `MQRS_ltm12p5`        | mapping-quality rank-sum |
| `QD < 2.0`              | `QD_lt2`              | quality by depth         |
| `SOR > 3.0`             | `SOR_gt3`             | strand odds ratio        |

These expressions live in the `GATK4_VARIANTFILTRATION` block of `conf/modules.config` and are **not** exposed as `--params`; change them there or override `ext.args` with `-c`.

### Stage 2 — VCFtools iterative filtering

Starting from the `PASS`-only `*.snps.gatkfilters` VCF, a chain of single-task VCFtools processes applies, in order:

1. `VCF_MINDP` — site mean-depth floor (`--min-meanDP`)
2. three rounds alternating `VCF_SITE_1..3` (site `--max-missing`) with `VCF_RMINDV_1..3` (drop individuals whose `F_MISS` exceeds the round's cutoff)
3. `VCFTOOLS_DEPTH_OUTLIER` — high-depth outlier cut at `mean + N·SD`
4. optional `VCF_FINAL` (a final `--max-missing`, only when `--vcf_max_missing_final` is set)
5. `VCF_RMFILTERED` (`--remove-filtered-all`) → the cleaned VCF
6. analysis subsets from the cleaned VCF: `VCF_BIALLELIC`, `VCF_MAC`, `VCF_THIN`, and `VCF_INGROUP` (drops samples flagged `outgroup=1`)

Compressed, tabix-indexed products (`.vcf.gz` + `.vcf.gz.tbi`) land in `results/vcftools/`. `results/vcftools/reports/` holds the per-round `*.imiss` missingness tables and `*.remove.txt` removal lists, plus `vcftools_filtering_summary.tsv` — a per-stage table of how many SNPs and individuals survive. That summary is also rendered in the MultiQC report as a **VCFtools filtering** section (a table plus SNP/individual trend bargraphs). Ingroup-only products (`*.ingroup.biallelic.vcf.gz`) are produced only when at least one sample is flagged `outgroup=1`.

The common VCFtools thresholds are exposed as params for convenience:

| Param                     | Default       | Description                                                                    |
| ------------------------- | ------------- | ------------------------------------------------------------------------------ |
| `--skip_vcftools`         | `false`       | Skip the whole VCFtools stage.                                                 |
| `--vcf_min_meandp`        | `5`           | Initial site mean-depth floor (`--min-meanDP`).                                |
| `--vcf_site_missing`      | `0.5,0.4,0.3` | Comma-separated site `--max-missing` schedule (one value per round).           |
| `--vcf_indiv_missing`     | `0.9,0.7,0.5` | Per-round individual `F_MISS` removal cutoffs (paired with the site schedule). |
| `--vcf_depth_sd`          | `2`           | `--max-meanDP` = mean + N·SD of site mean depth.                               |
| `--vcf_max_missing_final` | _(unset)_     | Optional final `--max-missing` after the depth cut.                            |
| `--vcf_mac`               | `3`           | `--mac` for the biallelic subset.                                              |
| `--vcf_thin`              | `10000`       | `--thin` for the biallelic subset.                                             |

### Tuning filters via `ext.args`

Every per-process flag set lives in a `withName:` block in [`conf/modules.config`](../conf/modules.config). The VCFtools thresholds are wired to the `--vcf_*` params through `ext.args` closures, for example:

```groovy
withName: VCF_SITE_1 { ext.args = { "--max-missing ${params.vcf_site_missing.tokenize(',')[0]}" } }
withName: VCF_MAC    { ext.args = { "--mac ${params.vcf_mac}" } }
```

There are two ways to change filtering behaviour:

- **Adjust an exposed threshold** — set the relevant param on the command line or in a `-params-file`:

  ```bash
  nextflow run tileylab/pgp ... --vcf_mac 5 --vcf_site_missing 0.6,0.5,0.4 --vcf_min_meandp 8
  ```

- **Change the actual flags** — for options that aren't param-exposed (the GATK filter expressions, or adding/swapping a VCFtools flag), override the process `ext.args` in a custom config passed with `-c`:

  ```groovy
  // my_filters.config
  process {
      withName: GATK4_VARIANTFILTRATION {
          ext.args = [
              '--filter-expression "QD < 3.0"  --filter-name "QD_lt3"',
              '--filter-expression "FS > 50.0" --filter-name "FS_gt50"',
              '--filter-expression "MQ < 40.0" --filter-name "MQ_lt40"'
          ].join(' ')
      }
      withName: VCF_THIN { ext.args = '--thin 5000' }
  }
  ```

  ```bash
  nextflow run tileylab/pgp -profile docker -c my_filters.config --input ... --outdir ...
  ```

  A `-c` override wins over the in-repo `conf/modules.config`, so you can re-tune filters without editing the pipeline.

> [!NOTE]
> The missingness loop is fixed at three rounds (`VCF_SITE_1..3` / `VCF_RMINDV_1..3`), so `--vcf_site_missing` and `--vcf_indiv_missing` must each carry three comma-separated values. Changing the _number_ of rounds is a `workflows/pgp.nf` edit, not a config change.

## Polyploid genotyping (EBG + updog)

The candidate biallelic SNP set (ingroup, before MAC filtering and thinning) is genotyped under two polyploid models — [EBG](https://github.com/pblischak/polyploid-genotyping) and [updog](https://github.com/dcgerard/updog) — and the two are compared. Each sample's `ploidy` comes from the samplesheet; samples are **grouped by ploidy** and EBG (`diseq`) and updog are run once per ploidy group, so mixed-ploidy cohorts are supported. A single per-locus sequencing-error vector (a `samtools mpileup` over the ingroup BAMs, converted to mean Phred error) is shared across groups and fed to EBG.

Disable the stage with `--no_polyploids` to run only the standard front end (mapping → GATK → hard filtering) for ordinary diploid genotyping.

| Param             | Default | Description                                                                                   |
| ----------------- | ------- | --------------------------------------------------------------------------------------------- |
| `--no_polyploids` | `false` | Skip EBG/updog genotyping (diploid front end only).                                           |
| `--ebg_model`     | `diseq` | EBG model: `diseq` (H-W disequilibrium) or `hwe`.                                             |
| `--ebg_iters`     | `1000`  | EBG ECM iterations.                                                                           |
| `--updog_model`   | `norm`  | updog prior model (see `updog::multidog`: `norm`, `hw`, `bb`, `s1`, `f1`, `flex`, `uniform`). |

Outputs land under `results/genotyping/`: `ebg/` and `updog/` hold the per-ploidy dosage calls, and `genotype_concordance.tsv` + `genotype_concordance_summary.tsv` report where the two methods agree and disagree. EBG reports ALT-allele dosage and updog reports reference-allele dosage; the comparison harmonises both to ALT dosage before scoring agreement.

#### Tuning EBG and updog via `ext.args`

The common options are param-exposed (`--ebg_model`, `--ebg_iters`, `--updog_model`). For anything else, both genotyping modules accept an `ext.args` pass-through — set it in a `-c` config exactly like the GATK/VCFtools steps:

```groovy
// my_genotyping.config
process {
    // Extra flags appended to the `ebg <model>` command (see `ebg diseq --help`)
    withName: EBG_DISEQ {
        ext.args = '--tol 1e-12 --stop 1e-6'
    }
    // Extra named arguments forwarded to updog::multidog (R name = value, comma-separated)
    withName: UPDOG {
        ext.args = 'seq = 0.01, update_bias = FALSE'
    }
}
```

```bash
nextflow run tileylab/pgp -profile docker -c my_genotyping.config --input ... --outdir ...
```

EBG's `ext.args` are inserted before `--prefix` (the `-n/-l/-p/-t/-a/-e` flags are managed by the module); updog's `ext.args` are parsed as R named arguments and passed to `multidog` via `do.call`, so any [`multidog`](https://dcgerard.github.io/updog/reference/multidog.html) argument (`seq`, `bias_init`, `prior_vec`, `update_bias`, …) can be set without editing the pipeline.

> [!NOTE]
> The genotyping images (`gptiley/ebg`, `gptiley/updog`) must contain `procps` (the `ps` command), which Nextflow uses to collect per-task metrics. The pinned images already do — if you swap in a different genotyping container, make sure it includes `procps`, or the process fails immediately with "Command 'ps' required by nextflow … cannot be found".

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
