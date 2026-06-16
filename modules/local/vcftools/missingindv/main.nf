process VCFTOOLS_MISSING_INDV {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/vcftools:0.1.16--pl5321hdcf5f25_11'
        : 'quay.io/biocontainers/vcftools:0.1.16--pl5321hdcf5f25_11'}"

    input:
    tuple val(meta), path(vcf)

    output:
    tuple val(meta), path("${prefix}.vcf"), emit: vcf
    path "${prefix}.imiss"                , emit: imiss
    path "${prefix}.remove.txt"           , emit: removed
    tuple val("${task.process}"), val('vcftools'), eval("vcftools --version 2>&1 | head -n1 | sed 's/VCFtools (//;s/)//'"), topic: versions, emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix        = task.ext.prefix ?: "${meta.id}"
    def threshold = task.ext.args ?: '1.0'
    """
    vcftools_missing_indv.sh ${vcf} ${prefix} ${threshold}
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.vcf ${prefix}.imiss ${prefix}.remove.txt
    """
}
