include { GFF_TSEBRA_SPFILTERFEATUREFROMKILLLIST                    } from '../../subworkflows/local/gff_tsebra_spfilterfeaturefromkilllist'
include { GFFCOMPARE as COMPARE_BRAKER_TO_LIFTOFF                   } from '../../modules/nf-core/gffcompare/main'
include { AGAT_SPFILTERFEATUREFROMKILLLIST                          } from '../../modules/nf-core/agat/spfilterfeaturefromkilllist/main'
include { GFFCOMPARE as VALIDATE_PURGING_BY_AGAT                    } from '../../modules/nf-core/gffcompare/main'
include { AGAT_SPMERGEANNOTATIONS as MERGE_BRAKER_LIFTOFF           } from '../../modules/nf-core/agat/spmergeannotations/main'

workflow PURGE_BRAKER_MODELS {
    take:
    braker_gff3                 // [ meta, gff3 ]
    braker_hints                // [ meta, gff ]
    liftoff_gff3                // [ meta, gff3 ]
    tsebra_config               // Channel: [ cfg ]
    allow_isoforms              // val(true|false)

    main:
    ch_versions                 = Channel.empty()

    // SUBWORKFLOW: GFF_TSEBRA_SPFILTERFEATUREFROMKILLLIST
    GFF_TSEBRA_SPFILTERFEATUREFROMKILLLIST(
        braker_gff3,
        braker_hints,
        tsebra_config,
        allow_isoforms,
        'braker'
    )

    ch_tsebra_killed_gff        = GFF_TSEBRA_SPFILTERFEATUREFROMKILLLIST.out.tsebra_killed_gff
    ch_versions                 = ch_versions.mix(GFF_TSEBRA_SPFILTERFEATUREFROMKILLLIST.out.versions)

    // MODULE: GFFCOMPARE as COMPARE_BRAKER_TO_LIFTOFF
    ch_comparison_inputs        = ch_tsebra_killed_gff
                                | join(liftoff_gff3)


    COMPARE_BRAKER_TO_LIFTOFF (
        ch_comparison_inputs.map { meta, braker, liftoff -> [ meta, braker ] },
        [ [], [], [] ],
        ch_comparison_inputs.map { meta, braker, liftoff -> [ meta, liftoff ] },
    )

    ch_tracking                 = COMPARE_BRAKER_TO_LIFTOFF.out.tracking
    ch_versions                 = ch_versions.mix(COMPARE_BRAKER_TO_LIFTOFF.out.versions.first())

    // COLLECTFILE: Transcript level kill list
    ch_kill_list                = ch_tracking
                                | map { meta, tracking ->

                                    def kept_lines = tracking.readLines()
                                        .findAll { line ->
                                            def cols = line.split('\t')

                                            ( cols[3] != 'u' ) && ( cols[3] != 'p' )
                                        }

                                    def tx_kill_list = kept_lines
                                        .collect { line ->
                                            def cols = line.split('\t')

                                            def matched = cols[4] =~ /q1:([^\|]+)\|([^\|]+)/

                                            matched[0][2].trim()
                                        }.join('\n')

                                    [ "${meta.id}.kill.list.txt" ] + tx_kill_list
                                }
                                | collectFile(newLine: true)
                                | map { file ->
                                    [ [ id: file.baseName.replace('.kill.list', '') ], file ]
                                }

    // MODULE: AGAT_SPFILTERFEATUREFROMKILLLIST
    ch_agat_kill_inputs         = ch_tsebra_killed_gff
                                | join(ch_kill_list)


    AGAT_SPFILTERFEATUREFROMKILLLIST(
        ch_agat_kill_inputs.map { meta, gff, kill -> [ meta, gff ] },
        ch_agat_kill_inputs.map { meta, gff, kill -> kill },
        [] // default config
    )

    ch_braker_purged_gff        = AGAT_SPFILTERFEATUREFROMKILLLIST.out.gff
    ch_versions                 = ch_versions.mix(AGAT_SPFILTERFEATUREFROMKILLLIST.out.versions.first())

    // Handle case where liftoff is not present
    ch_all_braker_gff           = ch_tsebra_killed_gff
                                | join(ch_braker_purged_gff, remainder:true)
                                | map { meta, tsebra_gff, purged_gff ->
                                    if ( purged_gff ) { return [ meta, purged_gff ] }
                                    if ( tsebra_gff ) { return [ meta, tsebra_gff ] }
                                }

    emit:
    braker_purged_gff           = ch_all_braker_gff     // [ meta, gff3 ]
    versions                    = ch_versions           // [ versions.yml ]
}
