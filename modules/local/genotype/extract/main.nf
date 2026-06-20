process GENOTYPE_EXTRACT {
    tag "ploidy ${ploidy}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container 'docker.io/python:3.12'

    input:
    tuple val(meta), path(vcf), val(ploidy)
    path ploidy_map

    output:
    tuple val(ploidy), path("p${ploidy}.tot.txt"), path("p${ploidy}.alt.txt")     , emit: ebg
    tuple val(ploidy), path("p${ploidy}.ref.txt"), path("p${ploidy}.size.txt")    , emit: updog
    tuple val(ploidy), path("p${ploidy}.samples.txt"), path("p${ploidy}.loci.txt"), emit: ids
    tuple val("${task.process}"), val('python'), eval("python3 --version | sed 's/Python //'"), topic: versions, emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    vcf_to_genotype_matrices.py \\
        --vcf ${vcf} \\
        --ploidy-map ${ploidy_map} \\
        --ploidy ${ploidy} \\
        --prefix p${ploidy}
    """

    stub:
    """
    printf '0 0\\n' > p${ploidy}.tot.txt
    printf '0 0\\n' > p${ploidy}.alt.txt
    printf 'snp\\ts1\\nchr1_1\\t0\\n' > p${ploidy}.ref.txt
    printf 'snp\\ts1\\nchr1_1\\t0\\n' > p${ploidy}.size.txt
    printf 's1\\n' > p${ploidy}.samples.txt
    printf 'chr1\\t1\\n' > p${ploidy}.loci.txt
    """
}
