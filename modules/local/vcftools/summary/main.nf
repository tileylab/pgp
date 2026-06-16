process VCFTOOLS_SUMMARY {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/vcftools:0.1.16--pl5321hdcf5f25_11'
        : 'quay.io/biocontainers/vcftools:0.1.16--pl5321hdcf5f25_11'}"

    input:
    tuple val(meta), val(labels), path(vcfs)

    output:
    path "vcftools_filtering_summary.tsv", emit: summary
    path "*_mqc.{tsv,yaml}"              , emit: mqc

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    vcftools_filter_summary.sh "${labels.join(',')}" ${vcfs}
    """

    stub:
    """
    printf 'stage\\tindividuals\\tsnps\\n' > vcftools_filtering_summary.tsv
    touch vcftools_filtering_table_mqc.tsv vcftools_filtering_snps_mqc.yaml vcftools_filtering_individuals_mqc.yaml
    """
}
