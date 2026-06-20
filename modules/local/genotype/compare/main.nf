process GENOTYPE_COMPARE {
    tag "${meta.id}"
    label 'process_single'

    // Dedicated minimal python image, independent of the EBG/updog containers
    // so a future change to either tool cannot perturb the comparison.
    conda "${moduleDir}/environment.yml"
    container 'docker.io/python:3.12'

    input:
    tuple val(meta), path(group_files, stageAs: 'inputs/*')

    output:
    path "genotype_concordance.tsv"        , emit: concordance
    path "genotype_concordance_summary.tsv", emit: summary
    tuple val("${task.process}"), val('python'), eval("python3 --version | sed 's/Python //'"), topic: versions, emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    # Per-group inputs (ebg_p*-genos.txt, p*.samples.txt, p*.loci.txt, updog_p*.updog.tsv)
    # are collected into inputs/; compare_genotypes.py discovers groups by ploidy.
    cp inputs/* .
    compare_genotypes.py
    """

    stub:
    """
    printf 'snp\\tind\\tploidy\\tebg_alt_dosage\\tupdog_alt_dosage\\tstatus\\n' > genotype_concordance.tsv
    printf 'metric\\tvalue\\noverall_concordance\\tNA\\n' > genotype_concordance_summary.tsv
    """
}
