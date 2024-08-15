#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include { validateParameters    } from 'plugin/nf-validation'

validateParameters()

include { GENEPAL               } from './workflows/genepal.nf'

workflow {
    PFR_GENEPAL()
}

workflow PFR_GENEPAL {
    GENEPAL()
}
