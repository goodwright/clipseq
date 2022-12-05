process BEDTOOLS_SHIFT {
    tag "$meta.id"
    label 'process_single'

    conda (params.enable_conda ? "bioconda::bedtools=2.30.0" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bedtools:2.30.0--hc088bd4_0' :
        'quay.io/biocontainers/bedtools:2.30.0--hc088bd4_0' }"

    input:
    tuple val(meta), path(bed)
    path(fai)

    output:
    tuple val(meta), path("*.bed"), emit: bed
    path  "versions.yml"          , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: '-s 0'
    def prefix = task.ext.prefix ?: "${meta.id}"
    def name = task.ext.suffix ? "${prefix}${task.ext.suffix}" : "${prefix}.shifted"
    """
    bedtools \\
        shift \\
        $args \\
        -i $bed \\
        -g $fai \\
        > ${name}.bed

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bedtools: \$(bedtools --version | sed -e "s/bedtools v//g")
    END_VERSIONS
    """
}
