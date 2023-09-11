process PARACLU_PARACLU {
    tag "$meta.id"
    label "process_low"

    conda "bioconda::paraclu=10"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/paraclu:10--h9a82719_1' :
        'biocontainers/paraclu:10--h9a82719_1' }"

    input:
    tuple val(meta), path(bed)
    val min_value

    output:
    tuple val(meta), path("*.sigxls.tsv"), emit: tsv
    path "versions.yml"                  , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args    = task.ext.args ? "$args" : "$min_value"
    def prefix  = task.ext.prefix ?: "${meta.id}"
    def VERSION = '10' // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    """
    paraclu $args $bed > ${prefix}.sigxls.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        paraclu: $VERSION
    END_VERSIONS
    """
}
