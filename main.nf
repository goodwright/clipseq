#!/usr/bin/env nextflow
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    goodwright/clipseq
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Github : https://github.com/goodwright/clipseq
----------------------------------------------------------------------------------------
*/

nextflow.enable.dsl = 2

include { goodwright_logo } from './modules/goodwright/util/logging/main'

workflow {
    log.info goodwright_logo()
}