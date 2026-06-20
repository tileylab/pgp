process EBG_DISEQ {
    tag "ploidy ${ploidy}"
    label 'process_low'

    // EBG is a custom binary with no conda package; container required.
    container 'docker.io/gptiley/ebg:main-55e543f'

    input:
    tuple val(ploidy), path(tot), path(alt), path(error)

    output:
    tuple val(ploidy), path("ebg_p${ploidy}-genos.txt"), emit: genos
    path "ebg_p${ploidy}-{PL,freqs}.txt"               , emit: reports
    tuple val("${task.process}"), val('ebg'), eval("/app/ebg --version 2>&1 | sed -n 's/.*ebg version \\([^ ]*\\).*/\\1/p'"), topic: versions, emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def iters = params.ebg_iters
    def model = params.ebg_model
    def args  = task.ext.args ?: ''
    """
    N=\$(wc -l < ${tot})
    L=\$(awk 'NR==1{print NF}' ${tot})
    /app/ebg ${model} \\
        -n \$N -l \$L -p ${ploidy} \\
        -t ${tot} -a ${alt} -e ${error} \\
        --iters ${iters} \\
        ${args} \\
        --prefix ebg_p${ploidy}
    """

    stub:
    """
    printf '0\\t0\\n' > ebg_p${ploidy}-genos.txt
    printf '0\\n'     > ebg_p${ploidy}-PL.txt
    printf '0\\n'     > ebg_p${ploidy}-freqs.txt
    """
}
