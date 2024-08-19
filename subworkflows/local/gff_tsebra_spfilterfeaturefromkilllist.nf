include { AGAT_CONVERTSPGFF2GTF                                     } from '../../modules/nf-core/agat/convertspgff2gtf/main'
include { TSEBRA                                                    } from '../../modules/nf-core/tsebra/main'
include { AGAT_CONVERTSPGXF2GXF                                     } from '../../modules/nf-core/agat/convertspgxf2gxf/main'
include { AGAT_SPFILTERFEATUREFROMKILLLIST as KILL_TSEBRA_ISOFORMS  } from '../../modules/nf-core/agat/spfilterfeaturefromkilllist/main'

workflow GFF_TSEBRA_SPFILTERFEATUREFROMKILLLIST {

    take:
    input_gff3                  // [ meta, gff3 ]
    braker_hints                // [ meta, gff ]
    tsebra_config               // Channel: [ cfg ]
    allow_isoforms              // val(true|false)
    val_prefix                  // val(String)

    main:
    ch_versions                 = Channel.empty()

    // MODULE: AGAT_CONVERTSPGFF2GTF
    AGAT_CONVERTSPGFF2GTF ( input_gff3 )

    ch_input_gtf                = AGAT_CONVERTSPGFF2GTF.out.output_gtf
    ch_versions                 = ch_versions.mix(AGAT_CONVERTSPGFF2GTF.out.versions.first())

    // COLLECTFILE: Prepare for TSEBRA
    ch_tsebra_input_gtf         = ch_input_gtf
                                | map { meta, gtf ->

                                    def lines = gtf.readLines()
                                        .collect { line ->
                                            if ( line.startsWith('#') ) { return line }

                                            def cols = line.split('\t')
                                            def feat = cols[2]

                                            if ( ! ( feat in [ 'gene', 'transcript', 'mRNA' ] ) ) { return line }

                                            def atts    = cols[8]
                                            def matches = atts =~ /ID ([^;]*)/
                                            def id      = matches[0][1]

                                            def feat_format = ( feat == 'mRNA' ) ? 'transcript' : feat

                                            return ( cols[0..1] + [ feat_format ] + cols[3..7] + [ id ] ).join('\t')
                                        }.join('\n')

                                    [ "${meta.id}.clean.gtf" ] + [ lines ]
                                }
                                | collectFile(newLine: true)
                                | map { file ->
                                    [ [ id: file.baseName.replace(".clean", "") ], file ]
                                }

    // MODULE: TSEBRA
    ch_tsebra_inputs            = ch_tsebra_input_gtf
                                | join(braker_hints)
                                | combine(tsebra_config)
    TSEBRA(
        ch_tsebra_inputs.map { meta, gtf, gff, cfg -> [ meta, [ gtf ] ] },
        ch_tsebra_inputs.map { meta, gtf, gff, cfg -> [ gff ] },
        [],
        ch_tsebra_inputs.map { meta, gtf, gff, cfg -> cfg }
    )

    ch_tsebra_gtf               = TSEBRA.out.tsebra_gtf
    ch_versions                 = ch_versions.mix(TSEBRA.out.versions.first())

    // COLLECTFILE: Format TSEBRA output
    ch_tsebra_formatted_gtf     = ch_tsebra_gtf
                                | map { meta, gtf ->

                                    def lines = gtf.readLines()
                                        .collect { line ->
                                            if ( line.startsWith('#') ) { return line }

                                            def cols    = line.split('\t')
                                            def program = cols[1]
                                            def feat    = cols[2]
                                            def atts    = cols[8]

                                            def atts_r  = ''
                                            // Remove attributes and use AGAT_CONVERTSPGXF2GXF
                                            // to create attributes based on sequential layout

                                            if ( feat != 'transcript' || program != 'Liftoff' ) {
                                                return ( cols[0..7] + [ atts_r ] ).join('\t')
                                            }

                                            def tx_id = atts.trim().replaceFirst('anno1.', '')
                                            def atts_g  = "liftoffID $tx_id"

                                            return ( cols[0..7] + [ atts_g ] ).join('\t')
                                        }.join('\n')

                                    [ "${meta.id}.gtf" ] + [ lines ]
                                }
                                | collectFile(newLine: true)
                                | map { file ->
                                    [ [ id: file.baseName ], file ]
                                }

    // MODULE: AGAT_CONVERTSPGXF2GXF
    AGAT_CONVERTSPGXF2GXF ( ch_tsebra_formatted_gtf )

    ch_tsebra_formatted_gff     = AGAT_CONVERTSPGXF2GXF.out.output_gff
    ch_versions                 = ch_versions.mix(AGAT_CONVERTSPGXF2GXF.out.versions.first())

    // COLLECTFILE: Format AGAT_CONVERTSPGXF2GXF output
    ch_tsebra_gff               = ch_tsebra_formatted_gff
                                | map { meta, gff ->

                                    def lines = gff.readLines()
                                        .collect { line ->
                                            if ( line.startsWith('#') ) { return line }

                                            def cols    = line.split('\t')
                                            def atts_r  = cols[8].replaceAll('-', '').replaceAll('agat', '')

                                            return ( cols[0..7] + [ atts_r ] ).join('\t')
                                        }.join('\n')

                                    [ "${meta.id}.${val_prefix}.gff3" ] + [ lines ]
                                }
                                | collectFile(newLine: true)
                                | map { file ->
                                    [ [ id: file.baseName.replace(".${val_prefix}", '') ], file ]
                                }

    // COLLECTFILE: Iso-form kill list if allow_isoforms=true
    ch_post_tsebra_kill_list    = allow_isoforms
                                ? Channel.empty()
                                : ch_tsebra_gff
                                | map { meta, gff ->
                                    def kill_list = gff.readLines()
                                        .findAll { line ->
                                            if ( line.startsWith('#') ) { return false }

                                            def cols    = line.split('\t')
                                            def feat    = cols[2]

                                            ( feat == 'mRNA' || feat == 'transcript' )
                                        }
                                        .collect { line ->
                                            def cols    = line.split('\t')
                                            def atts    = cols[8]
                                            def tx_id   = ( atts =~ /ID=([^;]*)/ )[0][1]
                                            def g_id    = ( atts =~ /Parent=([^;]*)/ )[0][1]

                                            [ g_id, tx_id ]
                                        }
                                        .groupBy { g_id, tx_id -> g_id }
                                        .findAll { key, value -> value.size() > 1 }
                                        .collect { key, value ->
                                            value.collect { it[1] }[1..-1]
                                        }
                                        .flatten()
                                        .join('\n')

                                    [ "${meta.id}.kill.list.txt" ] + [ kill_list ]
                                }
                                | collectFile(newLine: true)
                                | map { file ->
                                    [ [ id: file.baseName.replace('.kill.list', '') ], file ]
                                }

    // MODULE: AGAT_SPFILTERFEATUREFROMKILLLIST as KILL_TSEBRA_ISOFORMS
    ch_tsebra_kill_inputs       = ch_tsebra_gff
                                | join(ch_post_tsebra_kill_list)


    KILL_TSEBRA_ISOFORMS(
        ch_tsebra_kill_inputs.map { meta, gff, kill -> [ meta, gff ] },
        ch_tsebra_kill_inputs.map { meta, gff, kill -> kill },
        [] // default config
    )

    ch_tsebra_killed_gff        = ch_tsebra_gff
                                | join(KILL_TSEBRA_ISOFORMS.out.gff, remainder: true)
                                | map { meta, tsebra, killed ->
                                    if ( tsebra ) { [ meta, killed ?: tsebra ] }
                                }
    ch_versions                 = ch_versions.mix(KILL_TSEBRA_ISOFORMS.out.versions.first())

    emit:
    tsebra_killed_gff           = ch_tsebra_killed_gff  // [ val(meta), gff ]
    versions                    = ch_versions           // [ versions.yml ]
}
