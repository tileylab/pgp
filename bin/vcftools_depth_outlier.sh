#!/usr/bin/env bash
# Remove sites whose mean depth exceeds mean + k*SD (a high-depth outlier cut).
#
# Usage: vcftools_depth_outlier.sh <in.vcf> <out_prefix> <k_sd>
#
# Computes per-site mean depth (--site-mean-depth), derives the cutoff
# mean + k*SD, and applies --max-meanDP. With too little data to compute a
# cutoff the input is passed through unchanged. Writes <out_prefix>.vcf.
set -euo pipefail

in_vcf=$1
out_prefix=$2
k_sd=$3

vcftools --vcf "$in_vcf" --site-mean-depth --out "${out_prefix}_depth"

max_dp=$(awk -v k="$k_sd" '
    NR>1 && $3!="" { x+=$3; xx+=$3*$3; n++ }
    END {
        if (n>1)      { mu=x/n; var=(xx-n*mu*mu)/(n-1); if (var<0) var=0; printf "%.4f", mu+k*sqrt(var) }
        else if (n==1){ printf "%.4f", x }
    }' "${out_prefix}_depth.ldepth.mean")

if [ -n "$max_dp" ]; then
    vcftools --vcf "$in_vcf" --max-meanDP "$max_dp" --recode --recode-INFO-all --out "$out_prefix"
    mv "${out_prefix}.recode.vcf" "${out_prefix}.vcf"
else
    cp "$in_vcf" "${out_prefix}.vcf"
fi
