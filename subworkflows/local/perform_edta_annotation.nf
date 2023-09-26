nextflow.enable.dsl=2

include { SHORTEN_EDTA_IDS  } from '../../modules/local/edta/shorten_edta_ids'
include { EDTA              } from '../../modules/local/edta/edta'
include { RESTORE_EDTA_IDS  } from '../../modules/local/edta/restore_edta_ids'

// https://github.com/Plant-Food-Research-Open/assembly_qc
// GPL-3.0: https://github.com/Plant-Food-Research-Open/assembly_qc/blob/main/LICENSE
workflow PERFORM_EDTA_ANNOTATION {
    take:
        genome_fasta    // [meta, /path/to/genome/fasta]
    
    main:
        genome_fasta
        | SHORTEN_EDTA_IDS

        SHORTEN_EDTA_IDS
        .out
        .renamed_ids_fasta
        | EDTA

        RESTORE_EDTA_IDS(
            EDTA.out.te_anno_gff3,
            EDTA.out.intact_gff3.map { it[1] },
            EDTA.out.pass_list.map { it[1] },
            EDTA.out.out_file.map { it[1] },
            EDTA.out.te_lib_fasta.map { it[1] },
            SHORTEN_EDTA_IDS.out.renamed_ids_tsv.map { it[1] }
        )
    
    emit:
        te_anno_gff3    = RESTORE_EDTA_IDS.out.te_anno_gff3
        intact_gff3     = RESTORE_EDTA_IDS.out.intact_gff3
        pass_list       = RESTORE_EDTA_IDS.out.pass_list
        out_file        = RESTORE_EDTA_IDS.out.out_file
        te_lib_fasta    = RESTORE_EDTA_IDS.out.te_lib_fasta
        renamed_ids_tsv = RESTORE_EDTA_IDS.out.RESTORE_EDTA_IDS
}