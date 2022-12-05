//
// Formats crosslink file in BED format for input into paraclu then runs paraclu and paraclu-cut before reformatting the output
//

/*
* MODULES
*/
include { PARACLU_PARACLU                  } from '../../../modules/goodwright/paraclu/paraclu/main.nf'
include { PARACLU_CUT                      } from '../../../modules/goodwright/paraclu/cut/main.nf'
include { LINUX_COMMAND as PARACLU_PREPARE } from '../../../modules/goodwright/linux/command/main.nf'
include { LINUX_COMMAND as PARACLU_CONVERT } from '../../../modules/goodwright/linux/command/main.nf'

workflow PARACLU_ANALYSE {
    take:
    bed               // channel: [ val(meta), [ bed ] ]
    paraclu_min_value // val: integer

    main:
    ch_versions = Channel.empty()

    /*
    * MODULE: Prepare bedfile for input into paraclu using AWK
    */
    PARACLU_PREPARE (
        bed,
        [],
        false
    )
    ch_versions = ch_versions.mix(PARACLU_PREPARE.out.versions)

    /*
    * MODULE: Run paraclu (finds clusters in data attached to sequences)
    */
    PARACLU_PARACLU (
        PARACLU_PREPARE.out.file,
        paraclu_min_value
    )
    ch_versions = ch_versions.mix(PARACLU_PARACLU.out.versions)

    /*
    * MODULE: Run paraclu-cut (subset the output of paraclu)
    */
    PARACLU_CUT (
        PARACLU_PARACLU.out.tsv
    )
    ch_versions = ch_versions.mix(PARACLU_PARACLU.out.versions)

    /*
    * MODULE: Re-format parclu-cut output
    */
    PARACLU_CONVERT (
        PARACLU_CUT.out.tsv,
        [],
        false
    )

    emit:
    sigxls   = PARACLU_PARACLU.out.tsv  // channel: [ val(meta), [ tsv ] ]
    peaks    = PARACLU_CONVERT.out.file // channel: [ val(meta), [ tsv ] ]
    versions = ch_versions              // channel: [ versions.yml ]
}
