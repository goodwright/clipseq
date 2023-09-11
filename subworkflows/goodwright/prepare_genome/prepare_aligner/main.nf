//
// Uncompresses and prepares aligner indexes
//

include { UNTAR as UNTAR_BT1  } from '../../../../modules/nf-core/untar/main'
include { UNTAR as UNTAR_STAR } from '../../../../modules/nf-core/untar/main'
include { BOWTIE_BUILD        } from '../../../../modules/nf-core/bowtie/build/main'
include { STAR_GENOMEGENERATE } from '../../../../modules/nf-core/star/genomegenerate/main'

workflow PREPARE_ALINGER {
    take:
    aligners        // list: aligners, any of [bowtie2, star]
    fasta           // channel: [ fasta ]
    gtf             // channel: [ gtf ]
    bt1_index_path  // channel: [ folder/tar.gz ]
    star_index_path // channel: [ folder/tar.gz ]

    main:
    ch_versions = Channel.empty()

    // Check for meta data channel and add if not present
    ch_fasta_input = fasta
    if(!(fasta instanceof groovyx.gpars.dataflow.DataflowVariable ||
       fasta instanceof groovyx.gpars.dataflow.DataflowBroadcast)) {
        ch_fasta_input = Channel.of([[:], fasta])
    }

    // Check for meta data channel and add if not present
    ch_gtf_input = gtf
    if(!(gtf instanceof groovyx.gpars.dataflow.DataflowVariable || 
       gtf instanceof groovyx.gpars.dataflow.DataflowBroadcast)) {
        ch_gtf_input = Channel.of([[:], gtf])
    }

    /*
    * MODULES: Uncompress Bowtie2 index or generate if required
    */
    ch_bt1_index = Channel.empty()
    if ("bowtie" in aligners && bt1_index_path) {
        if (bt1_index_path.toString().endsWith(".tar.gz")) {
            ch_bt1_index = UNTAR_BT1 ( [ [:], bt1_index_path ] ).untar
            ch_versions  = ch_versions.mix(UNTAR_BT1.out.versions)
        } else {
            ch_bt1_index = Channel.of([[:], bt1_index_path])
        }
    }
    else if ("bowtie" in aligners) {
        ch_bt1_index = BOWTIE_BUILD ( ch_fasta_input.collect{it[1]} ).index.map{ [[:], it] }
        ch_versions  = ch_versions.mix(BOWTIE_BUILD.out.versions)
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
            ch_star_index = Channel.of([[:], star_index_path])
        }
    }
    else if ("star" in aligners) {
        ch_star_index = STAR_GENOMEGENERATE ( ch_fasta_input, ch_gtf_input ).index
        ch_versions   = ch_versions.mix(STAR_GENOMEGENERATE.out.versions)
    }

    emit:
    bt_index   = ch_bt1_index  // channel: [ val(meta), [ folder ] ]
    star_index = ch_star_index // channel: [ val(meta), [ folder ] ]
    versions   = ch_versions   // channel: [ versions.yml ]
}
