process GENOTYPE_ERROR {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container 'docker.io/python:3.12'

    input:
    tuple val(meta), path(pileup), path(positions)

    output:
    path "error.txt", emit: error
    tuple val("${task.process}"), val('python'), eval("python3 --version | sed 's/Python //'"), topic: versions, emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    per_locus_error.py --pileup ${pileup} --loci ${positions} --out error.txt
    """

    stub:
    """
    printf '0.01\\n' > error.txt
    """
}
