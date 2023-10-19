nextflow.enable.dsl=2

include { GUNZIP as GUNZIP_TARGET_ASSEMBLY      } from '../modules/nf-core/gunzip'
include { GUNZIP as GUNZIP_TE_LIBRARY           } from '../modules/nf-core/gunzip'
include { GUNZIP as GUNZIP_EXTERNAL_PROTEIN_SEQ } from '../modules/nf-core/gunzip'
include { FASTA_VALIDATE                        } from '../modules/local/fasta_validate'
include { REPEATMASKER                          } from '../modules/kherronism/repeatmasker'
include { STAR_GENOMEGENERATE                   } from '../modules/nf-core/star/genomegenerate'
include { CAT_FASTQ                             } from '../modules/nf-core/cat/fastq'
include { SORTMERNA                             } from '../modules/nf-core/sortmerna'
include { STAR_ALIGN                            } from '../modules/nf-core/star/align'
include { SAMTOOLS_CAT                          } from '../modules/nf-core/samtools/cat'
include { CAT_CAT as CAT_PROTEIN_SEQS           } from '../modules/nf-core/cat/cat'
include { BRAKER3                               } from '../modules/kherronism/braker3'

include { PERFORM_EDTA_ANNOTATION               } from '../subworkflows/local/perform_edta_annotation'
include { EXTRACT_SAMPLES                       } from '../subworkflows/local/extract_samples'
include { FASTQ_FASTQC_UMITOOLS_FASTP           } from '../subworkflows/nf-core/fastq_fastqc_umitools_fastp'

include { validateParams                        } from '../modules/local/validate_params'

validateParams(params)

// Additional validation
// Check rRNA databases for sortmerna
if (params.sample_prep.remove_ribo_rna) {
    ch_ribo_db = file(params.sample_prep.ribo_database_manifest, checkIfExists: true)
    if (ch_ribo_db.isEmpty()) {exit 1, "File provided with --ribo_database_manifest is empty: ${ch_ribo_db.getName()}!"}
}

workflow PAN_GENE {

    // Versions
    Channel.empty()
    | set { ch_versions }
    
    // GUNZIP: target_assemblies
    Channel.fromList(params.target_assemblies)
    | map { tag, filePath ->
        [[id:tag], file(filePath, checkIfExists: true)]
    }
    | branch { meta, file ->
        gz: "$file".endsWith(".gz")
        rest: !"$file".endsWith(".gz")
    }
    | set { ch_target_assemblies }

    GUNZIP_TARGET_ASSEMBLY(
        ch_target_assemblies.gz
    )
    .gunzip
    | mix(
        ch_target_assemblies.rest
    )
    | set { ch_gunzip_target_assemblies }

    ch_versions
    | mix(GUNZIP_TARGET_ASSEMBLY.out.versions.first())
    | set { ch_versions }

    // FASTA_VALIDATE
    FASTA_VALIDATE(ch_gunzip_target_assemblies)
    .valid_fasta
    | set { ch_validated_target_assemblies }

    ch_versions
    | mix(FASTA_VALIDATE.out.versions.first())
    | set { ch_versions }

    // GUNZIP: te_libraries
    Channel.fromList(params.te_libraries)
    | map { tag, filePath ->
        [[id:tag], file(filePath, checkIfExists: true)]
    }
    | branch { meta, file ->
        gz: "$file".endsWith(".gz")
        rest: !"$file".endsWith(".gz")
    }
    | set { ch_te_libraries }

    GUNZIP_TE_LIBRARY(
        ch_te_libraries.gz
    )
    .gunzip
    | mix(
        ch_te_libraries.rest
    )
    | set { ch_gunzip_te_libraries }

    ch_versions
    | mix(GUNZIP_TE_LIBRARY.out.versions.first())
    | set { ch_versions }

    // PERFORM_EDTA_ANNOTATION
    ch_validated_target_assemblies
    | join(
        ch_gunzip_te_libraries, remainder: true
    )
    | filter { meta, assembly, teLib ->
        teLib == null
    }
    | map {meta, assembly, teLib -> [meta, assembly]}
    | PERFORM_EDTA_ANNOTATION

    ch_versions
    | mix(PERFORM_EDTA_ANNOTATION.out.versions)
    | set { ch_versions }
    
    // REPEATMASKER
    ch_validated_target_assemblies
    | join(
        PERFORM_EDTA_ANNOTATION.out.te_lib_fasta.mix(ch_gunzip_te_libraries)
    )
    | set { ch_assemblies_n_te_libs }

    REPEATMASKER(
        ch_assemblies_n_te_libs.map {meta, assembly, teLib -> [meta, assembly]},
        ch_assemblies_n_te_libs.map {meta, assembly, teLib -> teLib},
    )

    ch_versions
    | mix(REPEATMASKER.out.versions.first())
    | set { ch_versions }

    // STAR_GENOMEGENERATE
    def star_ignore_sjdbgtf = true
    STAR_GENOMEGENERATE(
        REPEATMASKER.out.fasta_masked,
        REPEATMASKER.out.fasta_masked.map{meta, maskedFasta -> [meta, []]},
        star_ignore_sjdbgtf
    )
    .index
    | set { ch_assembly_index }

    ch_versions
    | mix(STAR_GENOMEGENERATE.out.versions.first())
    | set { ch_versions }

    // EXTRACT_SAMPLES
    // https://github.com/nf-core/rnaseq
    // MIT: https://github.com/nf-core/rnaseq/blob/master/LICENSE
    // Changes
    // Use meta.id as key for groupTuple as groupTuple does not work when there is a sublist in the key list
    ch_samplesheet_path = Channel.empty()
    if(params.samplesheet != null) {
        ch_samplesheet_path = Channel.fromPath(params.samplesheet)
    }
    
    EXTRACT_SAMPLES(
        ch_samplesheet_path,
        Channel.of(params.target_assemblies.collect{tag, fastaPath -> tag.strip()}.join(","))
    )
    .reads
    | map { meta, fastq ->
        new_id = meta.id - ~/_T\d+/
        [ new_id, meta + [id: new_id], fastq ]
    }
    | groupTuple()
    | branch { meta_id, meta, fastq ->
        single  : fastq.size() == 1
            return [ meta.first(), fastq.flatten() ]
        multiple: fastq.size() > 1
            return [ meta.first(), fastq.flatten() ]
    }
    | set { ch_fastq }

    ch_versions
    | mix(EXTRACT_SAMPLES.out.versions)
    | set { ch_versions }

    // CAT_FASTQ
    // https://github.com/nf-core/rnaseq
    // MIT: https://github.com/nf-core/rnaseq/blob/master/LICENSE
    CAT_FASTQ (
        ch_fastq.multiple
    )
    .reads
    | mix(ch_fastq.single)
    | set { ch_cat_fastq }
    
    ch_versions
    | mix(CAT_FASTQ.out.versions.first())
    | set { ch_versions }

    // FASTQ_FASTQC_UMITOOLS_FASTP
    // https://github.com/nf-core/rnaseq
    // MIT: https://github.com/nf-core/rnaseq/blob/master/LICENSE
    def with_umi            = false
    def skip_umi_extract    = true
    def umi_discard_read    = false
    FASTQ_FASTQC_UMITOOLS_FASTP (
        ch_cat_fastq,
        params.sample_prep.skip_fastqc,
        with_umi,
        skip_umi_extract,
        umi_discard_read,
        params.sample_prep.skip_fastp,
        [],
        params.sample_prep.save_trimmed,
        params.sample_prep.save_trimmed,
        params.sample_prep.min_trimmed_reads
    )
    .reads
    | set { ch_trim_reads }

    // SORTMERNA
    if (params.sample_prep.remove_ribo_rna) {
        Channel.from(ch_ribo_db.readLines())
        | map { row -> file(row, checkIfExists: true) }
        | collect
        | set { ch_sortmerna_fastas }

        SORTMERNA (
            ch_trim_reads,
            ch_sortmerna_fastas
        )
        .reads
        | set { ch_trim_reads }

        ch_versions
        | mix(SORTMERNA.out.versions.first())
        | set { ch_versions }
    }

    ch_trim_reads
    | flatMap { meta, reads ->
        def targetAssemblies = meta["target_assemblies"]

        readsByAssembly = []

        for(assembly in targetAssemblies) {
            readsByAssembly += [[[id: "${meta.id}.on.${assembly}", single_end: meta.single_end, target_assembly: assembly], reads]]
        }

        return readsByAssembly
    }
    | set { ch_trim_reads_by_assembly }

    ch_versions
    | mix(FASTQ_FASTQC_UMITOOLS_FASTP.out.versions)
    | set { ch_versions }

    // STAR_ALIGN
    ch_assembly_index
    | map { meta, index ->
        [meta.id, index]
    }
    | cross(
        ch_trim_reads_by_assembly.map{meta, reads -> [meta.target_assembly, meta, reads]}
    )
    | map { indexWithExt, readsWithExt ->
        def index = indexWithExt[1]

        def readsMeta = readsWithExt[1]
        def reads = readsWithExt[2]

        [
            readsMeta,
            reads,
            index
        ]
    }
    | set { ch_trim_reads_by_assembly_with_index }

    def seq_platform = false
    def seq_center = false
    STAR_ALIGN(
        ch_trim_reads_by_assembly_with_index.map{meta, reads, index -> [meta, reads]},
        ch_trim_reads_by_assembly_with_index.map{meta, reads, index -> [[id: meta.target_assembly], index]},
        ch_trim_reads_by_assembly_with_index.map{meta, reads, index -> [[id: meta.target_assembly], []]},
        star_ignore_sjdbgtf,
        seq_platform,
        seq_center
    )
    .bam_sorted
    .tap { ch_mapped_reads }
    .map { meta, bam ->
        [
            [id: meta.target_assembly],
            bam instanceof List ? bam.find {it =~ /Aligned/} : bam
        ]
    }
    | groupTuple
    | set { ch_mapped_reads_by_assembly }

    ch_versions
    | mix(STAR_ALIGN.out.versions.first())
    | set { ch_versions }

    // SAMTOOLS_CAT
    ch_mapped_reads_by_assembly
    | branch { meta, bamList ->
        bams: bamList.size() > 1
        bam: bamList.size() <= 1
    }
    | set { ch_samtools_cat_inputs_branches }

    SAMTOOLS_CAT(
        ch_samtools_cat_inputs_branches.bams
    )
    .bam
    | map { meta, bam ->
        [
            meta,
            [bam]
        ]
    }
    | mix(
        ch_samtools_cat_inputs_branches.bam
    )
    | set { ch_cat_bam_by_assembly }

    ch_versions
    | mix(SAMTOOLS_CAT.out.versions.first())
    | set { ch_versions }

    // GUNZIP: external_protein_seqs
    ch_external_protein_seqs = Channel.empty()
    if(params.external_protein_seqs != null) {
        ch_external_protein_seqs = Channel.fromList(params.external_protein_seqs)
    }
    
    ch_external_protein_seqs
    | map { filePath ->
        def fileHandle = file(filePath, checkIfExists: true)
        [[id:fileHandle.getSimpleName()], fileHandle]
    }
    | branch { meta, file ->
        gz: "$file".endsWith(".gz")
        rest: !"$file".endsWith(".gz")
    }
    | set { ch_external_protein_seqs_branch }

    GUNZIP_EXTERNAL_PROTEIN_SEQ(
        ch_external_protein_seqs_branch.gz
    )
    .gunzip
    | mix(
        ch_external_protein_seqs_branch.rest
    )
    | set { ch_gunzip_external_protein_seqs }

    ch_versions
    | mix(GUNZIP_EXTERNAL_PROTEIN_SEQ.out.versions.first())
    | set { ch_versions }

    // CAT_PROTEIN_SEQS
    ch_gunzip_external_protein_seqs
    | map{meta, filePath -> filePath}
    | collect
    | map{fileList -> [[id:"protein_seqs"], fileList]}
    | CAT_PROTEIN_SEQS
    
    CAT_PROTEIN_SEQS.out.file_out
    | set { ch_protein_seq }

    ch_versions
    | mix(CAT_PROTEIN_SEQS.out.versions)
    | set { ch_versions }

    // BRAKER3
    REPEATMASKER.out.fasta_masked
    | mix(ch_cat_bam_by_assembly)
    | groupTuple(size: 2, remainder: true)
    | map { meta, groupedItems ->
        def maskedFasta = groupedItems[0]

        if(groupedItems.size() == 2) {
            def bam = groupedItems[1]
            return [meta, maskedFasta, bam]
        } else {
            return [meta, maskedFasta, []]
        }
    }
    | set { ch_braker_inputs }
    
    if(params.external_protein_seqs) {
        ch_braker_inputs
        | combine(ch_protein_seq.map{meta, filePath -> filePath})
        | set { ch_braker_inputs }
    } else {
        ch_braker_inputs
        | map{meta, assembly, bams -> [meta, assembly, bams, []]}
        | set { ch_braker_inputs }
    }
    
    ch_fasta            = ch_braker_inputs.map{meta, assembly, bams, proteinSeq -> [meta, assembly]}
    ch_bam              = ch_braker_inputs.map{meta, assembly, bams, proteinSeq -> bams}
    ch_rnaseq_sets_dirs = ch_braker_inputs.map{meta, assembly, bams, proteinSeq -> []}
    ch_rnaseq_sets_ids  = ch_braker_inputs.map{meta, assembly, bams, proteinSeq -> []}
    ch_proteins         = ch_braker_inputs.map{meta, assembly, bams, proteinSeq -> proteinSeq}
    ch_hintsfile        = ch_braker_inputs.map{meta, assembly, bams, proteinSeq -> []}

    BRAKER3(
        ch_fasta,
        ch_bam,
        ch_rnaseq_sets_dirs,
        ch_rnaseq_sets_ids,
        ch_proteins,
        ch_hintsfile
    )

    ch_versions
    | mix(BRAKER3.out.versions.first())
    | set { ch_versions }
}