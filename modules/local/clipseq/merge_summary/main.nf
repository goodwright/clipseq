process MERGE_SUMMARY {
    tag "$gtf"
    label "process_single"

    conda "conda-forge::pandas=1.4.3"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pandas:1.4.3':
        'biocontainers/pandas:1.4.3' }"

    input:
    tuple val(meta), path(summary_type)
    tuple val(meta), path(summary_subtype)
    tuple val(meta), path(summary_gene)
    tuple val(meta), path(smrna_premapped_k1_cDNA)

    output:
    tuple val(meta), path("*summary_type_premapadjusted.tsv")   , emit: summary_type_adjusted
    tuple val(meta), path("*summary_subtype_premapadjusted.tsv"), emit: summary_subtype_adjusted
    tuple val(meta), path("*summary_gene_premapadjusted.tsv")   , emit: summary_gene_adjusted
    path "versions.yml"                                         , emit: versions

    when:
    task.ext.when == null || task.ext.when

    shell:
    process_name = task.process
    template 'merge_summary.py'
}

