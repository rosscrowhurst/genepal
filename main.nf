#!/usr/bin/env nextflow
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    plant-food-research-open/genepal
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Github : https://github.com/plant-food-research-open/genepal
----------------------------------------------------------------------------------------
*/

nextflow.enable.dsl = 2

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS / WORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { GENEPAL                   } from './workflows/genepal'
include { PIPELINE_INITIALISATION   } from './subworkflows/local/utils_nfcore_genepal_pipeline'
include { PIPELINE_COMPLETION       } from './subworkflows/local/utils_nfcore_genepal_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    NAMED WORKFLOWS FOR PIPELINE
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// WORKFLOW: Run main analysis pipeline depending on type of input
//
workflow PLANTFOODRESEARCHOPEN_GENEPAL {

    take:
    target_assembly
    tar_assm_str
    is_masked
    te_library
    braker_annotation
    braker_ex_asm_str
    rna_fq
    rna_bam
    rna_bam_by_assembly
    sortmerna_fastas
    ext_prot_fastas
    liftoff_fasta
    liftoff_gff
    tsebra_config
    orthofinder_pep

    main:
    //
    // WORKFLOW: Run pipeline
    //
    GENEPAL(
        target_assembly
        tar_assm_str
        is_masked
        te_library
        braker_annotation
        braker_ex_asm_str
        rna_fq
        rna_bam
        rna_bam_by_assembly
        sortmerna_fastas
        ext_prot_fastas
        liftoff_fasta
        liftoff_gff
        tsebra_config
        orthofinder_pep
    )

}
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow {

    main:

    //
    // SUBWORKFLOW: Run initialisation tasks
    //
    PIPELINE_INITIALISATION (
        params.version,
        params.help,
        params.validate_params,
        params.monochrome_logs,
        args,
        params.outdir,
        params.input
    )

    //
    // WORKFLOW: Run main workflow
    //
    PLANTFOODRESEARCHOPEN_GENEPAL(
        PIPELINE_INITIALISATION.out.target_assembly
        PIPELINE_INITIALISATION.out.tar_assm_str
        PIPELINE_INITIALISATION.out.is_masked
        PIPELINE_INITIALISATION.out.te_library
        PIPELINE_INITIALISATION.out.braker_annotation
        PIPELINE_INITIALISATION.out.braker_ex_asm_str
        PIPELINE_INITIALISATION.out.rna_fq
        PIPELINE_INITIALISATION.out.rna_bam
        PIPELINE_INITIALISATION.out.rna_bam_by_assembly
        PIPELINE_INITIALISATION.out.sortmerna_fastas
        PIPELINE_INITIALISATION.out.ext_prot_fastas
        PIPELINE_INITIALISATION.out.liftoff_fasta
        PIPELINE_INITIALISATION.out.liftoff_gff
        PIPELINE_INITIALISATION.out.tsebra_config
        PIPELINE_INITIALISATION.out.orthofinder_pep
    )

    //
    // SUBWORKFLOW: Run completion tasks
    //
    PIPELINE_COMPLETION (
        params.email,
        params.email_on_fail,
        params.plaintext_email,
        params.outdir,
        params.monochrome_logs,
        params.hook_url
    )
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
