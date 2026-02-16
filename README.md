# tileylab/pgp

<!--
[![Open in GitHub Codespaces](https://img.shields.io/badge/Open_In_GitHub_Codespaces-black?labelColor=grey&logo=github)](https://github.com/codespaces/new/tileylab/pgp)
[![GitHub Actions CI Status](https://github.com/tileylab/pgp/actions/workflows/nf-test.yml/badge.svg)](https://github.com/tileylab/pgp/actions/workflows/nf-test.yml)
[![GitHub Actions Linting Status](https://github.com/tileylab/pgp/actions/workflows/linting.yml/badge.svg)](https://github.com/tileylab/pgp/actions/workflows/linting.yml)[![Cite with Zenodo](http://img.shields.io/badge/DOI-10.5281/zenodo.XXXXXXX-1073c8?labelColor=000000)](https://doi.org/10.5281/zenodo.XXXXXXX)
[![nf-test](https://img.shields.io/badge/unit_tests-nf--test-337ab7.svg)](https://www.nf-test.com)

[![Nextflow](https://img.shields.io/badge/version-%E2%89%A525.04.0-green?style=flat&logo=nextflow&logoColor=white&color=%230DC09D&link=https%3A%2F%2Fnextflow.io)](https://www.nextflow.io/)
[![nf-core template version](https://img.shields.io/badge/nf--core_template-3.5.2-green?style=flat&logo=nfcore&logoColor=white&color=%2324B064&link=https%3A%2F%2Fnf-co.re)](https://github.com/nf-core/tools/releases/tag/3.5.2)
[![run with conda](http://img.shields.io/badge/run%20with-conda-3EB049?labelColor=000000&logo=anaconda)](https://docs.conda.io/en/latest/)
[![run with docker](https://img.shields.io/badge/run%20with-docker-0db7ed?labelColor=000000&logo=docker)](https://www.docker.com/)
[![run with singularity](https://img.shields.io/badge/run%20with-singularity-1d355c.svg?labelColor=000000)](https://sylabs.io/docs/)
[![Launch on Seqera Platform](https://img.shields.io/badge/Launch%20%F0%9F%9A%80-Seqera%20Platform-%234256e7)](https://cloud.seqera.io/launch?pipeline=https://github.com/tileylab/pgp)
-->
## Introduction

>[!WARNING]
> This project is not intended for use yet. I am migrating it over from a basic diploid genotype caller workflow and it will take time to iron out the issues.

**tileylab/pgp** is a referenced-based genotyping pipeline. A user provides a reference genome and samplesheet of paired-end Illumina reads. Reads are mapped to the reference with BWA and genotypes called with GATK. Some read quality statisitics are reported along the way and basic hard-filtering performed at the end, but users need to check that filtering is sufficient and appropriate. Candidate variants will eventually be genotyped again under a suite of polyploid models to find a high-confidence set.

### Some simplifying assumtions
1. Only one read pair per individual are allowed
   * Sometimes multiple libraries are sequenced for the same individual, which requires assigning different read groups then merging bams after alignment with BWA. This scenario is highly unlikely for small genomes and not accounted for here. Technical replicates of the same library should be able to be concatendated prior to genotyping without incident.
2. There are no intervals
   * Intervals can be defined for large genomes, these would be specific known regions of the reference genome with some biologial relevance, such as individual chromosomes. These intervals are used to *scatter-gather* genotyping operations such that each chromosome is an independent compute job. One could define mulitple intervals on a genome comprised of a single haploid chromosome, but the computational benefits deiminish with many small intervals as opposed to a few large ones. Namely, this creates a lot of extra IO and extra work on recombining the resulting *vcfs* into a single genome-wide *vcf*. Thus, intervals are ignored in favor of not spinning up many trivial compute requests.

### An important note
Some filtering of the multisample vcf is done based on the GATK best practices. This might not always be a good idea, so the vcfs at various stages are provided. Some additional hard filtering can be done based on depth and missing data cutoffs provided in the config, but this is a bit subjective and may require data exploration and judgement calls if losing many samples or important samples. Thus, the goal of this pipeline is to get the computational load of genotyping done and dealing with all of the GATK steps in a containerized envrionment. Some additional bespoke software for poylploid genotyping have also been containerized to ease implementation barriers. Some care in the evaluation of results on a case-by-case basis is needed though.


Development of **tileylab/pgp** borrowed heavily from [sarek](https://github.com/nf-core/sarek) to learn how to finesse various nf-core modules. Many of the GATK modules were copied over to `modules/local` and modified for narrower purposes.

## Usage
The input data is a comma-separated sample sheet that looks as follows:

`samplesheet.csv`:

```csv
sample,fastq_1,fastq_2,species,ploidy
SRR19908723,data/SRR19908723.R1.fastq.gz,data/SRR19908723.R2.fastq.gz,Vdarrowii,2
SRR19908722,data/SRR19908722.R1.fastq.gz,data/SRR19908722.R2.fastq.gz,Vdarrowii,2
SRR19908705,data/SRR19908705.R1.fastq.gz,data/SRR19908705.R2.fastq.gz,Vdarrowii,2
```
Relative paths are acceptable too provided the root is the nextflow launch directory. Here is a real example where the sample sheet and the reads are in `root/data`


A reference genome is expected. This is for reference alignment and genotyping. An interval file is also needed for GATK. The parameters can be defined from the command line with `reference` and `intervals`, respectively. Here is an example:
```bash
nextflow run pgp --input pgp_data/samplesheet.csv --outdir pgp_result --intervals pgp_data/ref/intervals.list --reference pgp_data/ref/ref.fa --phix pgp_data/db/bbduk/phix174_ill.ref.fa -profile docker
```
The reference fasta file is also has a default in the config and json as `data/ref/ref.fa` but this can be overwritten at the command line if need be. Similar to the samplesheet and fastq availability, one might be able to provide a relative path or need to provide a full path depending on the resource.

A script for preparing the *intervals* file is available at `bin/prepare_intervals.py`. It might be automated in the pipeline to create the intervals file from the reference genome, but for now it is simply an isolated step the user needs to do and think about.

Some notable options have been pre-configure in `conf\modules.config`:
 * Only properly paired reads with a mapping score of at least 20 are retained for genotyping with '-b -q 20 -f 2' -F 4' applied by samtools view
 * Some cpu allocations need to be changed here rather than editing the profiles in the module main.nf scripts because the linting will complain
 * HaplotypeCaller from GATK runs in GVCF mode. This would need to be changed to *BP_RESOLUTION* if genotype information for all sites, such as differentiating between invariant and data-deficient, is needed.

## Credits

tileylab/pgp was originally written by George P. Tiley.

## Contributions and Support

If you would like to contribute to this pipeline, please see the [contributing guidelines](.github/CONTRIBUTING.md).

## Citations

<!-- TODO nf-core: Add citation for pipeline after first release. Uncomment lines below and update Zenodo doi and badge at the top of this file. -->
<!-- If you use tileylab/pgp for your analysis, please cite it using the following doi: [10.5281/zenodo.XXXXXX](https://doi.org/10.5281/zenodo.XXXXXX) -->

<!-- TODO nf-core: Add bibliography of tools and data used in your pipeline -->

An extensive list of references for the tools used by the pipeline can be found in the [`CITATIONS.md`](CITATIONS.md) file.

This pipeline uses code and infrastructure developed and maintained by the [nf-core](https://nf-co.re) community, reused here under the [MIT license](https://github.com/nf-core/tools/blob/main/LICENSE).

> **The nf-core framework for community-curated bioinformatics pipelines.**
>
> Philip Ewels, Alexander Peltzer, Sven Fillinger, Harshil Patel, Johannes Alneberg, Andreas Wilm, Maxime Ulysse Garcia, Paolo Di Tommaso & Sven Nahnsen.
>
> _Nat Biotechnol._ 2020 Feb 13. doi: [10.1038/s41587-020-0439-x](https://dx.doi.org/10.1038/s41587-020-0439-x).
