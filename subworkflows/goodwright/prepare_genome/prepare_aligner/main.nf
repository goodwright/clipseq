//
// Uncompresses and prepares aligner indexes
//

include { UNTAR as UNTAR_BT2  } from '../../../../modules/nf-core/untar/main'
include { UNTAR as UNTAR_STAR } from '../../../../modules/nf-core/untar/main'
include { BOWTIE2_BUILD       } from '../../../../modules/nf-core/bowtie2/build/main'
include { STAR_GENOMEGENERATE } from '../../../../modules/nf-core/star/genomegenerate/main'

workflow PREPARE_ALINGER {
    take:
    aligners        // list: aligners, any of [bowtie2, star]
    fasta           // channel: [ fasta ]
    gtf             // channel: [ gtf ]
    bt2_index_path  // channel: [ folder/tar.gz ]
    star_index_path // channel: [ folder/tar.gz ]

    main:
    ch_versions = Channel.empty()

    /*
    * MODULES: Uncompress Bowtie2 index or generate if required
    */
    ch_bt2_index = Channel.empty()
    if ("bowtie2" in aligners && bt2_index_path) {
        if (bt2_index_path.toString().endsWith(".tar.gz")) {
            ch_bt2_index = UNTAR_BT2 ( [ [:], bt2_index_path ] ).untar
            ch_versions  = ch_versions.mix(UNTAR_BT2.out.versions)
        } else {
            ch_bt2_index = bt2_index_path
        }
    }
    else if ("bowtie2" in aligners) {
        ch_bt2_index = BOWTIE2_BUILD ( fasta ).index
        ch_versions  = ch_versions.mix(BOWTIE2_BUILD.out.versions)
    }

    /*
    * MODULES: Uncompress STAR index or generate if required
    */
    ch_star_index = Channel.empty()
    if ("star" in aligners && star_index_path) {
        if (star_index_path.toString().endsWith(".tar.gz")) {
            ch_star_index = UNTAR_STAR ( [ [:], star_index_path ] ).untar
            ch_versions  = ch_versions.mix(UNTAR_STAR.out.versions)
        } else {
            ch_star_index = star_index_path
        }
    }
    else if ("star" in aligners) {
        // Check for meta data channel and strip if present
        ch_star_input = fasta
        if(fasta instanceof groovyx.gpars.dataflow.DataflowVariable) {
            ch_star_input = fasta.map{ [it.last()] }
        }
                
        ch_star_index = STAR_GENOMEGENERATE ( ch_star_input, gtf ).index
        ch_versions   = ch_versions.mix(STAR_GENOMEGENERATE.out.versions)
    }

    emit:
    bt2_index  = ch_bt2_index  // channel: [ val(meta), [ folder ] ]
    star_index = ch_star_index // channel: [ val(meta), [ folder ] ]
    versions   = ch_versions   // channel: [ versions.yml ]
}
