nextflow.enable.dsl=2

include { GUNZIP            } from '../modules/nf-core/gunzip'
include { FASTA_VALIDATE    } from '../modules/local/fasta_validate'

include { validateParams    } from '../modules/local/validate_params'

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

    GUNZIP(
        ch_target_assemblies.gz
    )
    .gunzip
    | mix(
        ch_target_assemblies.rest
    )
    | set { ch_gunzip_target_assemblies }

    ch_versions.mix(GUNZIP.out.versions)
    | set { ch_versions }

    // FASTA_VALIDATE
    ch_gunzip_target_assemblies
    | FASTA_VALIDATE
    | set { ch_validated_target_assemblies }

    ch_versions.mix(FASTA_VALIDATE.out.versions)
    | set { ch_versions }
}