process BGZIP_TABIX {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/samtools:1.22.1--h96c455f_0'
        : 'biocontainers/samtools:1.22.1--h96c455f_0'}"

    input:
    tuple val(meta), path(vcfs)

    output:
    tuple val(meta), path("*.vcf.gz")    , emit: vcf
    tuple val(meta), path("*.vcf.gz.tbi"), emit: tbi
    tuple val("${task.process}"), val('tabix'), eval("tabix --version 2>&1 | head -n1 | sed 's/^tabix (htslib) //'"), topic: versions, emit: versions_tabix

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    for v in ${vcfs}; do
        bgzip -c "\$v" > "\${v}.gz"
        tabix -p vcf "\${v}.gz"
    done
    """

    stub:
    """
    for v in ${vcfs}; do
        touch "\${v}.gz" "\${v}.gz.tbi"
    done
    """
}
