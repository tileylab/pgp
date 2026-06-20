process SAMTOOLS_MPILEUP {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/samtools:1.22.1--h96c455f_0'
        : 'biocontainers/samtools:1.22.1--h96c455f_0'}"

    input:
    tuple val(meta), path(vcf)
    path bams
    path bais
    tuple val(meta2), path(fasta)
    path fai

    output:
    path "pileup.txt"   , emit: pileup
    path "positions.tsv", emit: positions
    tuple val("${task.process}"), val('samtools'), eval("samtools --version | head -n1 | sed 's/samtools //'"), topic: versions, emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    # Candidate loci (coordinate order matches the EBG matrices built from the same VCF)
    grep -v '^#' ${vcf} | cut -f1,2 > positions.tsv

    printf '%s\\n' ${bams} > bamlist.txt

    samtools mpileup \\
        --fasta-ref ${fasta} \\
        --positions positions.tsv \\
        --bam-list bamlist.txt \\
        --output pileup.txt
    """

    stub:
    """
    printf 'chr1\\t1\\n' > positions.tsv
    printf 'chr1\\t1\\tA\\t1\\t.\\tI\\n' > pileup.txt
    """
}
