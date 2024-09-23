//
// Subworkflow with functionality specific to the plant-food-research-open/genepal pipeline
//

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { UTILS_NFVALIDATION_PLUGIN } from '../../nf-core/utils_nfvalidation_plugin'
include { paramsSummaryMap          } from 'plugin/nf-validation'
include { fromSamplesheet           } from 'plugin/nf-validation'
include { UTILS_NEXTFLOW_PIPELINE   } from '../../nf-core/utils_nextflow_pipeline'
include { completionEmail           } from '../../nf-core/utils_nfcore_pipeline'
include { completionSummary         } from '../../nf-core/utils_nfcore_pipeline'
include { dashedLine                } from '../../nf-core/utils_nfcore_pipeline'
include { nfCoreLogo                } from '../../nf-core/utils_nfcore_pipeline'
include { imNotification            } from '../../nf-core/utils_nfcore_pipeline'
include { UTILS_NFCORE_PIPELINE     } from '../../nf-core/utils_nfcore_pipeline'
include { workflowCitation          } from '../../nf-core/utils_nfcore_pipeline'

/*
========================================================================================
    SUBWORKFLOW TO INITIALISE PIPELINE
========================================================================================
*/

workflow PIPELINE_INITIALISATION {

    take:
    version           // boolean: Display version and exit
    help              // boolean: Display help text
    validate_params   // boolean: Boolean whether to validate parameters against the schema at runtime
    monochrome_logs   // boolean: Do not use coloured log outputs
    nextflow_cli_args //   array: List of positional nextflow CLI args
    outdir            //  string: The output directory where the results will be saved
    input             //  string: Path to input samplesheet

    main:

    //
    // Print version and exit if required and dump pipeline parameters to JSON file
    //
    UTILS_NEXTFLOW_PIPELINE (
        version,
        true,
        outdir,
        workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1
    )

    //
    // Validate parameters and generate parameter summary to stdout
    //
    pre_help_text = nfCoreLogo(monochrome_logs)
    post_help_text = '\n' + workflowCitation() + '\n' + dashedLine(monochrome_logs)
    def String workflow_command = """nextflow run ${workflow.manifest.name} \\
    -profile <docker/singularity/.../institute> \\
    --input assemblysheet.csv \\
    --protein_evidence proteins.faa \\
    --outdir <OUTDIR>"""
    UTILS_NFVALIDATION_PLUGIN (
        help,
        workflow_command,
        pre_help_text,
        post_help_text,
        validate_params,
        "nextflow_schema.json"
    )

    //
    // Check config provided to the pipeline
    //
    UTILS_NFCORE_PIPELINE (
        nextflow_cli_args
    )

    //
    // Create input channels
    //
    ch_input                    = Channel.fromSamplesheet('input')

    ch_target_assembly          = ch_input
                                | map { it ->
                                    def tag         = it[0]
                                    def fasta       = it[1]

                                    def fasta_file  = file(fasta, checkIfExists: true)
                                    def is_zipped   = fasta.endsWith('.gz')
                                    def sz_thresh   = is_zipped ? 300_000 : 1_000_000
                                    def fasta_size  = fasta_file.size()

                                    if ( fasta_size < sz_thresh ) { // < 1 MB
                                        error "The assembly represented by tag '$tag' is only $fasta_size bytes. The minimum allowed size is 1 MB!"
                                    }

                                    [ [ id: tag ], fasta_file ]
                                }

    ch_tar_assm_str             = ch_input
                                | map { it ->
                                    def tag         = it[0].strip()

                                    tag
                                }
                                | collect
                                | map { it ->
                                    it.join(",")
                                }

    ch_is_masked                = ch_input
                                | map { it ->
                                    def tag         = it[0]
                                    def is_masked   = it[2]

                                    [ [ id: tag ], is_masked == "yes" ]
                                }

    ch_te_library               = ch_input
                                | map { it ->
                                    def tag         = it[0]
                                    def te_fasta    = it[3]

                                    if ( te_fasta ) {
                                        [ [ id:tag ], file(te_fasta, checkIfExists: true) ]
                                    }
                                }

    ch_braker_annotation        = ch_input
                                | map { it ->
                                    def tag         = it[0]
                                    def braker_gff3 = it[4]
                                    def hints_gff   = it[5]

                                    if ( braker_gff3 ) {
                                        [
                                            [ id: tag ],
                                            file(braker_gff3, checkIfExists: true),
                                            file(hints_gff, checkIfExists: true)
                                        ]
                                    }
                                }

    ch_braker_ex_asm_str        = ch_braker_annotation
                                | map { meta, braker_gff3, hints_gff -> meta.id }
                                | collect
                                | map { it.join(",") }
                                | ifEmpty( "" )

    ch_rna_branch               = ! params.rna_evidence
                                ? Channel.empty()
                                : Channel.fromSamplesheet('rna_evidence')
                                | map { meta, f1, f2 ->
                                    f2
                                    ? [ meta + [ single_end: false ], [ file(f1, checkIfExists:true), file(f2, checkIfExists:true) ] ]
                                    : [ meta + [ single_end: true ], [ file(f1, checkIfExists:true) ] ]
                                }
                                | map { meta, files ->
                                    [ meta + [ target_assemblies: meta.target_assemblies.split(';').sort() ], files ]
                                }
                                | branch { meta, files ->
                                    fq:  files.first().extension != 'bam'
                                    bam: files.first().extension == 'bam'
                                }

    ch_rna_fq                   = ! params.rna_evidence
                                ? Channel.empty()
                                : ch_rna_branch.fq
                                | map { meta, files -> [ meta.id, meta, files ] }
                                | groupTuple
                                | combine(ch_tar_assm_str)
                                | map { id, metas, files, tar_assm_str ->
                                    validateFastqMetadata(metas, files, tar_assm_str)
                                }

    ch_rna_bam                  = ! params.rna_evidence
                                ? Channel.empty()
                                : ch_rna_branch.bam
                                | map { meta, files -> [ meta.id, meta, files ] }
                                | groupTuple
                                | combine(ch_tar_assm_str)
                                | flatMap { id, metas, files, tar_assm_str ->
                                    validateBamMetadata(metas, files, tar_assm_str)
                                }

    // Check if each sample for a given assembly has either bam or fastq files
    ch_rna_bam
    | flatMap { meta, bams ->
        meta.target_assemblies.collect { [ [ meta.id, it ], 'bam' ] }
    }
    | join(
        ch_rna_fq
        | flatMap { meta, fqs ->
            meta.target_assemblies.collect { [ [ meta.id, it ], 'fq' ] }
        }
    )
    | map { combination, bam, fq ->
        error "Sample ${combination[0]} for assembly ${combination[1]} can not have both fastq and bam files"
    }

    ch_rna_bam_by_assembly      = ch_rna_bam
                                | map { meta, bams -> [ [ id: meta.target_assemblies.first() ], bams ] }
                                | groupTuple
                                | map { meta, bams -> [ meta, bams.flatten() ] }

    ch_ribo_db                  = params.remove_ribo_rna
                                ? file(params.ribo_database_manifest, checkIfExists: true)
                                : null

    ch_sortmerna_fastas         = ch_ribo_db
                                ? Channel.from(ch_ribo_db ? ch_ribo_db.readLines() : null)
                                | map { row -> file(row, checkIfExists: true) }
                                | collect
                                : Channel.empty()

    ch_ext_prot_fastas          = ( params.protein_evidence.endsWith('txt')
                                    ? Channel.fromPath(params.protein_evidence)
                                    | splitText
                                    : Channel.fromPath(params.protein_evidence)
                                )
                                | map { file_path ->

                                    def file_handle = ( file_path instanceof String )
                                        ? file(file_path.strip(), checkIfExists: true)
                                        : file_path

                                    [ [ id: idFromFileName( file_handle.baseName ) ], file_handle ]
                                }


    ch_liftoff_mm               = ! params.liftoff_annotations
                                ? Channel.empty()
                                : Channel.fromSamplesheet('liftoff_annotations')
                                | multiMap { fasta, gff ->
                                    def fastaFile = file(fasta, checkIfExists:true)

                                    fasta: [ [ id: idFromFileName( fastaFile.baseName ) ], fastaFile ]
                                    gff: [ [ id: idFromFileName( fastaFile.baseName ) ], file(gff, checkIfExists:true) ]
                                }

    ch_liftoff_fasta            = params.liftoff_annotations
                                ? ch_liftoff_mm.fasta
                                : Channel.empty()

    ch_liftoff_gff              = params.liftoff_annotations
                                ? ch_liftoff_mm.gff
                                : Channel.empty()

    ch_tsebra_config            = Channel.of ( file("${projectDir}/assets/tsebra-template.cfg", checkIfExists: true) )
                                | map { cfg ->
                                    def param_intron_support = params.enforce_full_intron_support ? '1.0' : '0.0'

                                    def param_e1 = params.allow_isoforms ? '0.1'    : '0.0'
                                    def param_e2 = params.allow_isoforms ? '0.5'    : '0.0'
                                    def param_e3 = params.allow_isoforms ? '0.05'   : '0.0'
                                    def param_e4 = params.allow_isoforms ? '0.2'    : '0.0'

                                    [
                                        'tsebra-config.cfg',
                                        cfg
                                        .text
                                        .replace('PARAM_INTRON_SUPPORT', param_intron_support)
                                        .replace('PARAM_E1', param_e1)
                                        .replace('PARAM_E2', param_e2)
                                        .replace('PARAM_E3', param_e3)
                                        .replace('PARAM_E4', param_e4)
                                    ]
                                }
                                | collectFile


    ch_orthofinder_pep          = ! params.orthofinder_annotations
                                ? Channel.empty()
                                : Channel.fromSamplesheet('orthofinder_annotations')
                                | map { tag, fasta ->
                                    [ [ id: tag ], file(fasta, checkIfExists:true)  ]
                                }

    emit:
    target_assembly             = ch_target_assembly
    tar_assm_str                = ch_tar_assm_str
    is_masked                   = ch_is_masked
    te_library                  = ch_te_library
    braker_annotation           = ch_braker_annotation
    braker_ex_asm_str           = ch_braker_ex_asm_str
    rna_fq                      = ch_rna_fq
    rna_bam                     = ch_rna_bam
    rna_bam_by_assembly         = ch_rna_bam_by_assembly
    sortmerna_fastas            = ch_sortmerna_fastas
    ext_prot_fastas             = ch_ext_prot_fastas
    liftoff_fasta               = ch_liftoff_fasta
    liftoff_gff                 = ch_liftoff_gff
    tsebra_config               = ch_tsebra_config
    orthofinder_pep             = ch_orthofinder_pep
}

/*
========================================================================================
    SUBWORKFLOW FOR PIPELINE COMPLETION
========================================================================================
*/

workflow PIPELINE_COMPLETION {

    take:
    email           //  string: email address
    email_on_fail   //  string: email address sent on pipeline failure
    plaintext_email // boolean: Send plain-text email instead of HTML
    outdir          //    path: Path to output directory where results will be published
    monochrome_logs // boolean: Disable ANSI colour codes in log output
    hook_url        //  string: hook URL for notifications

    main:

    summary_params = paramsSummaryMap(workflow, parameters_schema: "nextflow_schema.json")

    //
    // Completion email and summary
    //
    workflow.onComplete {
        if (email || email_on_fail) {
            completionEmail(summary_params, email, email_on_fail, plaintext_email, outdir, monochrome_logs)
        }

        completionSummary(monochrome_logs)

        if (hook_url) {
            imNotification(summary_params, hook_url)
        }
    }

    workflow.onError {
        log.error "Pipeline failed. Please refer to troubleshooting docs: https://nf-co.re/docs/usage/troubleshooting"
    }
}

/*
========================================================================================
    FUNCTIONS
========================================================================================
*/

//
// Additional validation
//
def idFromFileName(fileName) {

    def trial = ( fileName
        ).replaceFirst(
            /\.f(ast)?q$/, ''
        ).replaceFirst(
            /\.f(asta|sa|a|as|aa|na)?$/, ''
        ).replaceFirst(
            /\.gff(3)?$/, ''
        ).replaceFirst(
            /\.gz$/, ''
        )

    if ( trial == fileName ) { return fileName }

    return idFromFileName ( trial )
}

def validateFastqMetadata(metas, fqs, permAssString) {
    def permAssList = permAssString.split(",")

    // Check if each listed assembly is permissible
    metas.each { meta ->
        if ( meta.target_assemblies.any { !permAssList.contains( it ) } ) {
            error "Sample ${meta.id} targets ${meta.target_assemblies} which are not in $permAssList"
        }
    }

    // Check if multiple runs of a sample have the same target assemblies
    if ( metas.collect { meta -> meta.target_assemblies }.unique().size() > 1 ) {
        error "Multiple runs of sample ${metas.first().id} must target same assemblies"
    }

    // Check if multiple runs of a sample have the same endedness
    if ( metas.collect { meta -> meta.single_end }.unique().size() > 1 ) {
        error "Multiple runs of sample ${metas.first().id} must have same endedness"
    }

    [ metas.first(), fqs ]
}


def validateBamMetadata(metas, bams, permAssString) {
    def permAssList = permAssString.split(",")

    // Check if each listed assembly is permissible
    metas.each { meta ->
        if ( meta.target_assemblies.any { !permAssList.contains( it ) } ) {
            error "Sample ${meta.id} targets ${meta.target_assemblies} which are not in $permAssList"
        }
    }

    // Check that when the first file is bam then the second file is absent
    bams.findAll { files ->
        files.first().extension == 'bam' && files.size() != 1
    }
    .each { error "Sample ${metas.first().id} contains both bam and fastq pairs. When a bam file is provided as file_1, a fastq for file_2 is not permitted" }

    // Check that a bam file only targets a single assembly
    bams.eachWithIndex { files, index ->
        if ( files.first().extension == 'bam' && metas[index].target_assemblies.size() > 1 ) {
            error "BAM file for sample ${metas.first().id} can only target one assembly: ${metas[index].target_assemblies}"
        }
    }

    metas.every { it.target_assemblies == metas.first().target_assemblies }
    ? [ [ metas.first(), bams.flatten() ] ]
    : metas.withIndex().collect { meta, index -> [ meta, bams[index].flatten() ] }
}
