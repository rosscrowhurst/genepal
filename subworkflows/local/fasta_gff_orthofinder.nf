include { FILE_GUNZIP as FASTA_GUNZIP   } from '../../subworkflows/local/file_gunzip'
include { FILE_GUNZIP as GFF_GUNZIP     } from '../../subworkflows/local/file_gunzip'
include { GFFREAD as EXTRACT_PROTEINS   } from '../../modules/nf-core/gffread/main'
include { ORTHOFINDER                   } from '../../modules/nf-core/orthofinder/main'

workflow FASTA_GFF_ORTHOFINDER {
    take:
    ch_pep_fasta                // [ meta, fasta ]
    ch_fasta                    // [ meta, fasta ]
    ch_gff                      // [ meta, gff ]

    main:
    ch_versions                 = Channel.empty()

    // SUBWORKFLOW: FILE_GUNZIP as FASTA_GUNZIP
    FASTA_GUNZIP ( ch_fasta )

    ch_fasta_unzipped           = FASTA_GUNZIP.out.gunzip
    ch_versions                 = ch_versions.mix(FASTA_GUNZIP.out.versions)

    // SUBWORKFLOW: FILE_GUNZIP as GFF_GUNZIP
    GFF_GUNZIP ( ch_gff )

    ch_gff_unzipped             = GFF_GUNZIP.out.gunzip
    ch_versions                 = ch_versions.mix(GFF_GUNZIP.out.versions)

    // MODULE: GFFREAD as EXTRACT_PROTEINS
    ch_extraction_inputs        = ch_gff_unzipped
                                | join(ch_fasta_unzipped)

    EXTRACT_PROTEINS(
        ch_extraction_inputs.map { meta, gff, fasta -> [ meta, gff ] },
        ch_extraction_inputs.map { meta, gff, fasta -> fasta }
    )

    ch_versions                 = ch_versions.mix(EXTRACT_PROTEINS.out.versions.first())

    // MODULE: ORTHOFINDER
    ch_orthofinder_peps         = EXTRACT_PROTEINS.out.gffread_fasta
                                | map { meta, fasta -> fasta }
                                | mix(
                                    ch_pep_fasta.map { meta, fasta -> fasta }
                                )
                                | collect
                                | filter { it.size() > 1 }

    ORTHOFINDER ( ch_orthofinder_peps.map { fastas -> [ [ id: 'pangene' ], fastas ] } )

    ch_versions                 = ch_versions.mix(ORTHOFINDER.out.versions)

    emit:
    versions                    = ch_versions           // [ versions.yml ]
}
