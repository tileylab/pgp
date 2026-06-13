process VCFTOOLS_FILTER {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/vcftools:0.1.16--pl5321hdcf5f25_11'
        : 'quay.io/biocontainers/vcftools:0.1.16--pl5321hdcf5f25_11'}"

    input:
    tuple val(meta), path(vcf), path(tbi)
    path outgroups

    output:
    tuple val(meta), path("${prefix}*.vcf")          , emit: products
    path "vcftools_filtering_summary.tsv"            , emit: summary
    path "*_mqc.{tsv,yaml}"                           , emit: mqc
    path "*.imiss"                                    , emit: imiss
    path "removed_individuals.tsv"                    , emit: removed
    tuple val("${task.process}"), val('vcftools'), eval("vcftools --version 2>&1 | head -n1 | sed 's/VCFtools (//;s/)//'"), topic: versions, emit: versions_vcftools

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix              = task.ext.prefix ?: "${meta.id}"
    def min_meandp      = params.vcf_min_meandp
    def site_list       = params.vcf_site_missing
    def indiv_list      = params.vcf_indiv_missing
    def depth_sd        = params.vcf_depth_sd
    def max_missing_fin = params.vcf_max_missing_final
    def mac             = params.vcf_mac
    def thin            = params.vcf_thin
    // Ingroup-only products are built only when an outgroup list was supplied
    def ingroup_block   = outgroups ? """
    if [ -s ${outgroups} ]; then
        vcftools --vcf "\$CLEAN" --remove ${outgroups} --min-alleles 2 --max-alleles 2 --recode --recode-INFO-all --out ${prefix}.ingroup.biallelic
        mv ${prefix}.ingroup.biallelic.recode.vcf ${prefix}.ingroup.biallelic.vcf
        record ingroup_biallelic ${prefix}.ingroup.biallelic.vcf
    fi
    """ : ''
    """
    # Provenance + per-stage summary files (always present)
    printf 'sample\\tround\\tcriterion\\n' > removed_individuals.tsv
    printf 'stage\\tindividuals\\tsnps\\n' > vcftools_filtering_summary.tsv

    # Record individuals (#CHROM columns - 9) and SNPs (non-# lines) for a stage VCF.
    # Stages get a 2-digit index prefix so order is preserved in the MultiQC table/plots.
    STAGE_N=0
    record() {  # \$1=label  \$2=vcf  \$3=optional 'gz'
        STAGE_N=\$((STAGE_N+1)); idx=\$(printf '%02d' \$STAGE_N)
        if [ "\${3:-}" = gz ]; then reader="gzip -cd"; else reader="cat"; fi
        ind=\$(\$reader "\$2" | awk '/^#CHROM/{print NF-9; exit}')
        snps=\$(\$reader "\$2" | grep -vc '^#' || true)
        printf '%s\\t%s\\t%s\\n' "\${idx}_\$1" "\${ind:-0}" "\${snps:-0}" >> vcftools_filtering_summary.tsv
    }

    record input ${vcf} gz

    # 1) Seed: drop sites below the mean-depth floor
    vcftools --gzvcf ${vcf} --min-meanDP ${min_meandp} --recode --recode-INFO-all --out step_dp
    CUR=step_dp.recode.vcf
    record min_meanDP "\$CUR"

    # 2) Iterative alternating site / individual missingness removal
    IFS=',' read -ra SITES <<< "${site_list}" || true
    IFS=',' read -ra INDIV <<< "${indiv_list}" || true
    n=\${#SITES[@]}
    for ((i=0; i<n; i++)); do
        s=\${SITES[\$i]}
        m=\${INDIV[\$i]}

        # site filter
        vcftools --vcf "\$CUR" --max-missing "\$s" --recode --recode-INFO-all --out site_\${i}
        CUR=site_\${i}.recode.vcf
        record "site_max_missing_\${s}" "\$CUR"

        # per-individual missingness -> list of individuals over the cutoff
        vcftools --vcf "\$CUR" --missing-indv --out round_\${i}
        awk -v T="\$m" 'NR>1 && (\$5+0) > (T+0) {print \$1}' round_\${i}.imiss > remove_\${i}.txt

        if [ -s remove_\${i}.txt ]; then
            awk -v r="\$i" -v T="\$m" '{print \$1"\\tsite_indiv_round_"r"\\tF_MISS>"T}' remove_\${i}.txt >> removed_individuals.tsv
            vcftools --vcf "\$CUR" --remove remove_\${i}.txt --recode --recode-INFO-all --out indiv_\${i}
            CUR=indiv_\${i}.recode.vcf
        fi
        record "indiv_F_MISS_\${m}" "\$CUR"
    done

    # 3) Depth-outlier cut: max-meanDP = mean + k*SD of site mean depth
    vcftools --vcf "\$CUR" --site-mean-depth --out sitedepth
    MAXDP=\$(awk -v k=${depth_sd} 'NR>1 && \$3!="" {x+=\$3; xx+=\$3*\$3; nn++} END{ if(nn>1){ mu=x/nn; var=(xx-nn*mu*mu)/(nn-1); if(var<0)var=0; printf "%.4f", mu+k*sqrt(var) } else if(nn==1){ printf "%.4f", x } }' sitedepth.ldepth.mean)
    if [ -n "\$MAXDP" ]; then
        vcftools --vcf "\$CUR" --max-meanDP "\$MAXDP" --recode --recode-INFO-all --out maxdepth
        CUR=maxdepth.recode.vcf
    fi
    record max_meanDP "\$CUR"

    # 4) Optional final site-missingness pass
    FINAL_MM="${max_missing_fin ?: ''}"
    if [ -n "\$FINAL_MM" ]; then
        vcftools --vcf "\$CUR" --max-missing "\$FINAL_MM" --recode --recode-INFO-all --out finalmm
        CUR=finalmm.recode.vcf
        record final_max_missing "\$CUR"
    fi

    # 5) Drop any sites still carrying a FILTER flag -> cleaned VCF
    vcftools --vcf "\$CUR" --remove-filtered-all --recode --recode-INFO-all --out ${prefix}.filtered
    mv ${prefix}.filtered.recode.vcf ${prefix}.filtered.vcf
    CLEAN=${prefix}.filtered.vcf
    record filtered "\$CLEAN"

    # 6) Analysis-ready subsets from the cleaned VCF
    vcftools --vcf "\$CLEAN" --min-alleles 2 --max-alleles 2 --recode --recode-INFO-all --out ${prefix}.biallelic
    mv ${prefix}.biallelic.recode.vcf ${prefix}.biallelic.vcf
    record biallelic ${prefix}.biallelic.vcf

    vcftools --vcf ${prefix}.biallelic.vcf --mac ${mac} --recode --recode-INFO-all --out ${prefix}.biallelic.mac${mac}
    mv ${prefix}.biallelic.mac${mac}.recode.vcf ${prefix}.biallelic.mac${mac}.vcf
    record biallelic_mac${mac} ${prefix}.biallelic.mac${mac}.vcf

    vcftools --vcf ${prefix}.biallelic.mac${mac}.vcf --thin ${thin} --recode --recode-INFO-all --out ${prefix}.biallelic.mac${mac}.thin${thin}
    mv ${prefix}.biallelic.mac${mac}.thin${thin}.recode.vcf ${prefix}.biallelic.mac${mac}.thin${thin}.vcf
    record biallelic_mac${mac}_thin${thin} ${prefix}.biallelic.mac${mac}.thin${thin}.vcf
    ${ingroup_block}

    # 7) MultiQC custom content: summary table + SNP/individual trend bargraphs
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
        printf 'Stage\\tIndividuals\\tSNPs\\n'
        tail -n +2 vcftools_filtering_summary.tsv
    } > vcftools_filtering_table_mqc.tsv

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
        awk -F'\\t' 'NR>1 {printf "  \\"%s\\":\\n    SNPs: %s\\n", \$1, \$3}' vcftools_filtering_summary.tsv
    } > vcftools_filtering_snps_mqc.yaml

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
        awk -F'\\t' 'NR>1 {printf "  \\"%s\\":\\n    Individuals: %s\\n", \$1, \$2}' vcftools_filtering_summary.tsv
    } > vcftools_filtering_individuals_mqc.yaml
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    printf 'sample\\tround\\tcriterion\\n' > removed_individuals.tsv
    printf 'stage\\tindividuals\\tsnps\\n' > vcftools_filtering_summary.tsv
    touch ${prefix}.filtered.vcf
    touch ${prefix}.biallelic.vcf
    touch ${prefix}.biallelic.mac${params.vcf_mac}.vcf
    touch ${prefix}.biallelic.mac${params.vcf_mac}.thin${params.vcf_thin}.vcf
    touch round_0.imiss
    touch vcftools_filtering_table_mqc.tsv
    touch vcftools_filtering_snps_mqc.yaml
    touch vcftools_filtering_individuals_mqc.yaml
    """
}
