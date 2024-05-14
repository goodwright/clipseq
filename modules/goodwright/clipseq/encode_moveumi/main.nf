process ENCODE_MOVEUMI {
    label "process_single"

    conda "bioconda::biopython=1.70"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
    'https://depot.galaxyproject.org/singularity/biopython:1.81 ' :
    'biocontainers/biopython:1.81' }"

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("*.fastq.gz"), emit: reads
    path  "versions.yml"               , emit: versions

    when:
    task.ext.when == null || task.ext.when

    shell:
    def args     = task.ext.args ?: ''
    prefix       = task.ext.prefix ?: "${meta.id}"
    process_name = task.process
    template 'encode_moveumi.py'
}
