nextflow.enable.dsl=2

include { validateParams } from '../modules/local/validate_params'

validateParams(params)

workflow PAN_GENE {
    Channel.fromList(params.target_assemblies)
}