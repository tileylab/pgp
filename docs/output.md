# tileylab/pgp: Output

## Introduction

This document describes the output produced by the pipeline. The directories listed below are created in the `--outdir` results directory after the pipeline finishes; all paths are relative to that top-level directory. VCFs are bgzip-compressed and tabix-indexed (`.vcf.gz` + `.vcf.gz.tbi`).

## Pipeline overview

The pipeline processes data through the following steps:

- [Read QC and preprocessing](#read-qc-and-preprocessing) — FastQC, fastp trimming, BBDuk PhiX decontamination
- [Alignment](#alignment) — BWA-MEM mapping, filtering, duplicate marking
- [Variant calling](#variant-calling) — GATK HaplotypeCaller + joint genotyping
- [Hard-filtered SNP sets](#hard-filtered-snp-sets) — GATK best-practice annotation filters
- [VCFtools filtering](#vcftools-filtering) — iterative missingness/depth filtering + analysis subsets
- [MultiQC](#multiqc) — aggregate QC report
- [Pipeline information](#pipeline-information) — run metrics and software versions

See [usage.md](usage.md#hard-filtering-strategy) for the full hard-filtering strategy and how to tune the thresholds.

### Read QC and preprocessing

<details markdown="1">
<summary>Output files</summary>

- `fastqc/` — `*_fastqc.html` / `*_fastqc.zip`: per-sample raw-read quality metrics.
- `fastp/` — `*.fastp.html` / `*.fastp.json` / `*.fastp.log`: adapter and quality trimming reports.
- `bbmap/` — `*.bbduk.log`: PhiX (sequencing spike-in) decontamination with BBDuk against `assets/phix174_ill.ref.fa`.

</details>

[FastQC](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/) reports per-base quality, GC content, adapter contamination, and overrepresented sequences. [fastp](https://github.com/OpenGene/fastp) trims adapters and low-quality bases; [BBDuk](https://jgi.doe.gov/data-and-tools/software-tools/bbtools/) removes PhiX reads before alignment.

### Alignment

<details markdown="1">
<summary>Output files</summary>

- `bwa/` — `*.bam`: reads mapped to the reference with BWA-MEM (read groups set from the sample id).
- `samtools/` — Q20, properly-paired BAMs (`samtools view -b -q 20 -f 2 -F 4`) plus `*.bai` indexes and `*.stats`/`*.flagstat`/`*.idxstats` summaries.
- `picard/` — `*.markdup.bam` + `*.MarkDuplicates.metrics.txt` (PCR duplicate marking; skipped when `--skip_markduplicates` is set).

</details>

### Variant calling

<details markdown="1">
<summary>Output files</summary>

- `gatk4/`
  - `<sample>.vcf.gz` — per-sample gVCFs from HaplotypeCaller (GVCF mode).
  - `joint_genotyping/` — the GenomicsDB workspace used for joint genotyping.
  - `joint_genotyping.vcf.gz` — the raw joint-genotyped multi-sample VCF.
  - `joint_genotyping.snps.annotations.vcf.gz` — SNPs with the GATK `FILTER` column populated by `VariantFiltration`.
  - `reference.dict` — sequence dictionary for the reference.

</details>

[GATK4](https://gatk.broadinstitute.org/) calls variants per sample in GVCF mode, then joint-genotypes the cohort. Because the target genomes are small, genotyping is run without interval scatter-gather.

### Hard-filtered SNP sets

<details markdown="1">
<summary>Output files</summary>

- `select/`
  - `joint_genotyping.snps.vcf.gz` — SNP sites only (`SELECT_ALL`).
  - `joint_genotyping.snps.gatkfilters.vcf.gz` — SNPs passing the GATK best-practice hard filters (`SELECT_PASS`); this is the input to the VCFtools stage.
  - `joint_genotyping.biallelic.gatkfilters.vcf.gz` — the PASS set restricted to biallelic SNPs (`SELECT_BIPASS`).

</details>

The GATK filter expressions (FS, MQ, ReadPosRankSum, MQRankSum, QD, SOR) are listed and explained in [usage.md](usage.md#stage-1--gatk-best-practice-hard-filters). The intermediate VCFs at each stage are kept so you can start from a less- or more-filtered set if the defaults are too aggressive.

### VCFtools filtering

<details markdown="1">
<summary>Output files</summary>

- `vcftools/` — final products of the iterative VCFtools stage:
  - `*.filtered.vcf.gz` — the cleaned VCF after the missingness loop, depth-outlier cut, and `--remove-filtered-all`.
  - `*.biallelic.vcf.gz`, `*.biallelic.mac<N>.vcf.gz`, `*.biallelic.mac<N>.thin<N>.vcf.gz` — analysis-ready biallelic / MAC-filtered / thinned subsets.
  - `*.ingroup.biallelic.vcf.gz` — biallelic set with outgroup samples removed (only when at least one sample is flagged `outgroup=1`).
- `vcftools/reports/`
  - `indiv_<n>.imiss` — per-individual missingness for each removal round.
  - `indiv_<n>.remove.txt` — individuals dropped in each round.
  - `vcftools_filtering_summary.tsv` — SNPs and individuals retained at every filtering stage.

</details>

The per-stage `vcftools_filtering_summary.tsv` is also rendered in the MultiQC report (a **VCFtools filtering** table plus SNP/individual trend bargraphs). The thresholds are tunable — see [usage.md](usage.md#stage-2--vcftools-iterative-filtering).

### Polyploid genotyping

<details markdown="1">
<summary>Output files</summary>

- `genotyping/ebg/` — `ebg_p<ploidy>-genos.txt` (ALT-allele dosage calls), `ebg_p<ploidy>-PL.txt`, `ebg_p<ploidy>-freqs.txt` per ploidy group.
- `genotyping/updog/` — `updog_p<ploidy>.updog.tsv` (per SNP × individual reference-allele dosage) per ploidy group.
- `genotyping/`
  - `genotype_concordance.tsv` — per SNP × individual, both calls (harmonised to ALT dosage) and an agreement flag.
  - `genotype_concordance_summary.tsv` — overall, per-ploidy and per-dosage-class agreement rates plus the most discordant SNPs.

</details>

The candidate biallelic SNP set is genotyped with [EBG](https://github.com/pblischak/polyploid-genotyping) (`diseq`) and [updog](https://github.com/dcgerard/updog), grouped by samplesheet ploidy, and the two are compared. Skip the stage with `--no_polyploids`. See [usage.md](usage.md#polyploid-genotyping-ebg--updog).

### MultiQC

<details markdown="1">
<summary>Output files</summary>

- `multiqc/`
  - `multiqc_report.html` — standalone HTML report.
  - `multiqc_data/` — parsed statistics behind the report.
  - `multiqc_plots/` — static images of the report plots.

</details>

[MultiQC](http://multiqc.info) aggregates QC across FastQC, fastp, BBDuk, samtools, Picard, and the VCFtools filtering summary into a single report, and records the software versions used.

### Pipeline information

<details markdown="1">
<summary>Output files</summary>

- `pipeline_info/`
  - Nextflow execution reports: `execution_report_*.html`, `execution_timeline_*.html`, `execution_trace_*.txt`, `pipeline_dag_*.html`.
  - `pgp_software_mqc_versions.yml` — software versions for every process.
  - `params_*.json` — the parameters used for the run.

</details>

[Nextflow](https://www.nextflow.io/docs/latest/tracing.html) generates the execution reports, which help troubleshoot failures and review run times and resource usage.
