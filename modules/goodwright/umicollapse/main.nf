process UMICOLLAPSE {
    tag "$meta.id"
    label "process_medium"

    container 'docker.io/elly1502/umicollapse:latest'

    input:
    tuple val(meta), path(bam), path(bai)

    output:
    tuple val(meta), path("*.bam"), emit: bam
    tuple val(meta), path("*.log"), emit: log
    path  "versions.yml"          , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    java -jar /UMICollapse/umicollapse.jar \\
        bam \\
        -i $bam \\
        -o ${prefix}.bam \\
        $args

    mv .command.log ${prefix}_UMICollapse.log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        umicollapse: NA
    END_VERSIONS
    """
}
