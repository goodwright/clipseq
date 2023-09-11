process ICOUNT_SEGMENT {
    tag "$gtf"
    label "process_single"

    // iCount-Mini segement did not work for Chris, but works in Charlotte's tests?
    // CC - I've put it back in because having different icount versions for segment & downstream functions
    // can cause weird issues sometimes.
    /* [RuntimeError] generator raised StopIteration
    File "/usr/local/lib/python3.9/site-packages/iCount/cli.py", line 448, in main
        result_object = func(**args)

    File "/usr/local/lib/python3.9/site-packages/iCount/genomes/segment.py", line 1015, in get_segments
        for gene_content in _get_gene_content(annotation, chromosomes, report_progress):
    */
    conda "bioconda::icount-mini=2.0.3"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/icount-mini:2.0.3--pyh5e36f6f_0' :
        'biocontainers/icount-mini:2.0.3--pyh5e36f6f_0' }"

    input:
    tuple val(meta), path(gtf)
    path fai

    output:
    tuple val(meta), path("*_seg.gtf")                  ,  emit: gtf
    tuple val(meta), path("*_regions.gtf.gz")           ,  emit: regions
    path "versions.yml"                                 ,  emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${gtf.simpleName}_seg"
    def regions_prefix = task.ext.regions_prefix ?: "${gtf.simpleName}"
    """
    iCount-Mini segment \\
        $gtf \\
        ${prefix}.gtf \\
        $fai

    mv regions.gtf.gz ${regions_prefix}_regions.gtf.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        iCount-Mini: \$(iCount-Mini -v)
    END_VERSIONS
    """
}
