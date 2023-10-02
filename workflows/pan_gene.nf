nextflow.enable.dsl=2

include { GUNZIP as GUNZIP_TARGET_ASSEMBLY  } from '../modules/nf-core/gunzip'
include { GUNZIP as GUNZIP_TE_LIBRARY       } from '../modules/nf-core/gunzip'
include { FASTA_VALIDATE                    } from '../modules/local/fasta_validate'
include { REPEATMASKER                      } from '../modules/kherronism/repeatmasker'
include { CAT_FASTQ                         } from '../modules/nf-core/cat/fastq'
include { BRAKER3                           } from '../modules/kherronism/braker3'

include { PERFORM_EDTA_ANNOTATION           } from '../subworkflows/local/perform_edta_annotation'
include { EXTRACT_SAMPLES                   } from '../subworkflows/local/extract_samples'
include { FASTQ_FASTQC_UMITOOLS_FASTP       } from '../subworkflows/nf-core/fastq_fastqc_umitools_fastp'

include { validateParams                    } from '../modules/local/validate_params'

validateParams(params)

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
    | mix(GUNZIP_TARGET_ASSEMBLY.out.versions)
    | set { ch_versions }

    // FASTA_VALIDATE
    FASTA_VALIDATE(ch_gunzip_target_assemblies)
    .valid_fasta
    | set { ch_validated_target_assemblies }

    ch_versions
    | mix(FASTA_VALIDATE.out.versions)
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
    | mix(GUNZIP_TE_LIBRARY.out.versions)
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
        ch_assemblies_n_te_libs.map {meta, assembly, te_lib -> [meta, assembly]},
        ch_assemblies_n_te_libs.map {meta, assembly, te_lib -> te_lib},
    )

    ch_versions
    | mix(REPEATMASKER.out.versions.first())
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

    ch_versions
    | mix(FASTQ_FASTQC_UMITOOLS_FASTP.out.versions)
    | set { ch_versions }

    // BRAKER3
    ch_bam = REPEATMASKER.out.fasta_masked.map {return []}
    ch_rnaseq_sets_dirs = REPEATMASKER.out.fasta_masked.map {return []}
    ch_rnaseq_sets_ids = REPEATMASKER.out.fasta_masked.map {return []}
    ch_proteins = REPEATMASKER.out.fasta_masked.map {return []}
    ch_hintsfile = REPEATMASKER.out.fasta_masked.map {return []}

    BRAKER3(
        REPEATMASKER.out.fasta_masked,
        ch_bam,
        ch_rnaseq_sets_dirs,
        ch_rnaseq_sets_ids,
        ch_proteins,
        ch_hintsfile
    )
}