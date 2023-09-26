nextflow.enable.dsl=2

include { GUNZIP            } from '../modules/nf-core/gunzip'

include { validateParams    } from '../modules/local/validate_params'

validateParams(params)

workflow PAN_GENE {
    
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
    | view
}