#!/usr/bin/env bash
# Summarise how many individuals and SNPs survive each filtering stage and emit
# a summary table plus MultiQC custom-content files.
#
# Usage: vcftools_filter_summary.sh "<label1,label2,...>" <vcf1> <vcf2> ...
#
# Labels and VCFs are paired positionally and must be in stage order. Writes:
#   vcftools_filtering_summary.tsv
#   vcftools_filtering_table_mqc.tsv
#   vcftools_filtering_snps_mqc.yaml
#   vcftools_filtering_individuals_mqc.yaml
set -euo pipefail

IFS=',' read -ra labels <<< "$1"
shift
vcfs=("$@")

# Read a VCF whether plain or gzipped.
read_vcf() {
    case "$1" in
        *.gz) gzip -cd "$1" ;;
        *)    cat "$1" ;;
    esac
}

printf 'stage\tindividuals\tsnps\n' > vcftools_filtering_summary.tsv
for i in "${!labels[@]}"; do
    vcf=${vcfs[$i]}
    # Single full-pass count: reading to EOF avoids SIGPIPE on the reader (an
    # early `awk ... exit` would kill `cat`/`gzip -cd` and trip `pipefail`).
    counts=$(read_vcf "$vcf" | awk '
        /^#CHROM/ { ind = NF - 9 }
        !/^#/     { snps++ }
        END       { printf "%d\t%d", ind+0, snps+0 }
    ')
    printf '%s\t%s\n' "${labels[$i]}" "$counts" >> vcftools_filtering_summary.tsv
done

# MultiQC: summary table
{
    cat <<'HEADER'
# id: "vcftools_filtering"
# section_name: "VCFtools filtering"
# description: "SNP sites and individuals retained at each stage of iterative VCFtools filtering."
# plot_type: "table"
# pconfig:
#     id: "vcftools_filtering_table"
#     namespace: "VCFtools filtering"
HEADER
    printf 'Stage\tIndividuals\tSNPs\n'
    tail -n +2 vcftools_filtering_summary.tsv
} > vcftools_filtering_table_mqc.tsv

# MultiQC: SNPs-retained bargraph
{
    cat <<'HEADER'
id: "vcftools_filtering_snps"
section_name: "VCFtools filtering: SNPs retained"
description: "SNP sites retained at each stage of iterative VCFtools filtering."
plot_type: "bargraph"
pconfig:
  id: "vcftools_filtering_snps_plot"
  title: "VCFtools filtering: SNPs retained per stage"
  ylab: "SNPs"
  cpswitch: false
data:
HEADER
    awk -F'\t' 'NR>1 {printf "  \"%s\":\n    SNPs: %s\n", $1, $3}' vcftools_filtering_summary.tsv
} > vcftools_filtering_snps_mqc.yaml

# MultiQC: individuals-retained bargraph
{
    cat <<'HEADER'
id: "vcftools_filtering_individuals"
section_name: "VCFtools filtering: individuals retained"
description: "Individuals retained at each stage of iterative VCFtools filtering."
plot_type: "bargraph"
pconfig:
  id: "vcftools_filtering_individuals_plot"
  title: "VCFtools filtering: individuals retained per stage"
  ylab: "Individuals"
  cpswitch: false
data:
HEADER
    awk -F'\t' 'NR>1 {printf "  \"%s\":\n    Individuals: %s\n", $1, $2}' vcftools_filtering_summary.tsv
} > vcftools_filtering_individuals_mqc.yaml
