process CLIPSEQ_FIND_LONGEST_TRANSCRIPT {
    tag "$gtf"
    label "process_single"

    conda     (params.enable_conda ? "conda-forge::python=3.10.4" : null)
    container "quay.io/biocontainers/python:3.10.4"

    input:
    tuple val(meta), path(gtf)

    output:
    tuple val(meta), path("*.txt")                 ,emit: longest_transcript
    tuple val(meta), path("*.fai")                 ,emit: longest_transcript_fai
    path  "versions.yml"                           ,emit: versions

    when:
    task.ext.when == null || task.ext.when

    shell:
    process_name = task.process
    output       = task.ext.output ?: "longest_transcript"
    template 'find_longest_transcript.py'
}
