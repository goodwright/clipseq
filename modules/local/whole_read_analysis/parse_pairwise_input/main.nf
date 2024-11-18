process PARSE_PAIRWISE_INPUT {
    tag "$samplesheet"
    label 'process_single'

    conda "conda-forge::python=3.10.4"
    container "docker.io/python:3.10.4"

    input:
    path samplesheet

    output:
    path '*.valid.csv'        , emit: csv
    path  "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    shell:
    process_name = task.process
    output       = task.ext.output ?: 'pairwise_samplesheet.valid.csv'
    template 'pairwise_samplesheet_check.py'
}
