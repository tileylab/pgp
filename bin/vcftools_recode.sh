#!/usr/bin/env bash
# Run a single VCFtools filtering operation and emit a plain VCF.
#
# Usage: vcftools_recode.sh <in.vcf[.gz]> <out_prefix> [vcftools args...]
#
# Detects a gzipped input automatically. Writes <out_prefix>.vcf.
set -euo pipefail

in_vcf=$1
out_prefix=$2
shift 2

case "$in_vcf" in
    *.gz) in_flag="--gzvcf" ;;
    *)    in_flag="--vcf" ;;
esac

vcftools "$in_flag" "$in_vcf" "$@" --recode --recode-INFO-all --out "$out_prefix"
mv "${out_prefix}.recode.vcf" "${out_prefix}.vcf"
