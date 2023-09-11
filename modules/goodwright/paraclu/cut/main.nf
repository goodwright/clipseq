process PARACLU_CUT {
    tag "$meta.id"
    label "process_single"

    conda "bioconda::paraclu=10"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/paraclu:10--h9a82719_1' :
        'biocontainers/paraclu:10--h9a82719_1' }"

    input:
    tuple val(meta), path(tsv)

    output:
    tuple val(meta), path("*.peaks.tsv"), emit: tsv
    path "versions.yml"                 , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args    = task.ext.args ?: ''
    def prefix  = task.ext.prefix ?: "${meta.id}"
    def VERSION = '10' // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    """
    paraclu-cut $args $tsv > ${prefix}.peaks.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        paraclu: $VERSION
    END_VERSIONS
    """
}
