#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include { PAN_GENE } from './workflows/pan_gene.nf'

workflow {
    PAN_GENE()
}