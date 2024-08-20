include { GUNZIP            } from '../../modules/nf-core/gunzip'

workflow FILE_GUNZIP {
    take:
    ch_input                // channel [ meta, archive ]

    main:

    ch_versions             = Channel.empty()

    ch_input_branch         = ch_input
                            | branch { meta, archive ->
                                gz: "$archive".endsWith('.gz')
                                rest: ! "$archive".endsWith('.gz')
                            }

    // MODULE: GUNZIP
    GUNZIP ( ch_input_branch.gz )

    ch_versions             = ch_versions.mix(GUNZIP.out.versions.first())

    emit:
    versions                = ch_versions
    gunzip                  = GUNZIP.out.gunzip.mix( ch_input_branch.rest )
}
