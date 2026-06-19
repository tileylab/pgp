#!/usr/bin/env bash
# Drop individuals whose genotype missingness exceeds a threshold.
#
# Usage: vcftools_missing_indv.sh <in.vcf> <out_prefix> <fmiss_threshold>
#
# Computes per-individual missingness (--missing-indv), lists individuals with
# F_MISS above the threshold, and removes them. If nobody exceeds the cutoff the
# input is passed through unchanged. Writes <out_prefix>.vcf, <out_prefix>.imiss,
# and <out_prefix>.remove.txt.
set -euo pipefail

in_vcf=$1
out_prefix=$2
threshold=$3

vcftools --vcf "$in_vcf" --missing-indv --out "$out_prefix"
awk -v T="$threshold" 'NR>1 && ($5+0) > (T+0) {print $1}' "${out_prefix}.imiss" > "${out_prefix}.remove.txt"

if [ -s "${out_prefix}.remove.txt" ]; then
    vcftools --vcf "$in_vcf" --remove "${out_prefix}.remove.txt" --recode --recode-INFO-all --out "$out_prefix"
    mv "${out_prefix}.recode.vcf" "${out_prefix}.vcf"
else
    cp "$in_vcf" "${out_prefix}.vcf"
fi
