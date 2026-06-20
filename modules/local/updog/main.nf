process UPDOG {
    tag "ploidy ${ploidy}"
    label 'process_medium'

    // updog runs in the maintainer's R container; no conda env used here.
    container 'docker.io/gptiley/updog:main-37c55ef'

    input:
    tuple val(ploidy), path(ref), path(size)

    output:
    tuple val(ploidy), path("updog_p${ploidy}.updog.tsv"), emit: genos
    tuple val("${task.process}"), val('updog'), eval("Rscript -e 'cat(as.character(packageVersion(\"updog\")))'"), topic: versions, emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def model = params.updog_model
    def args  = task.ext.args ?: ''
    """
    run_updog.R ${ref} ${size} ${ploidy} ${model} updog_p${ploidy} ${task.cpus} "${args}"
    """

    stub:
    """
    printf 'snp\\tind\\tgeno\\nchr1_1\\ts1\\t0\\n' > updog_p${ploidy}.updog.tsv
    """
}
