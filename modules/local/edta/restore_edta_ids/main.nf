nextflow.enable.dsl=2

// https://github.com/Plant-Food-Research-Open/assembly_qc
// GPL-3.0: https://github.com/Plant-Food-Research-Open/assembly_qc/blob/main/LICENSE
process RESTORE_EDTA_IDS {
    tag "${meta}"
    label "process_single"

    input:
        tuple val(meta), path(te_anno_gff3)
        path(intact_gff3)
        path(pass_list)
        path(out_file)
        path(te_lib_fa)
        path(renamed_ids_tsv)
    
    output:
        tuple val(meta), path("${meta}.EDTA.TEanno.gff3"),              emit: te_anno_gff3
        tuple val(meta), path("${meta}.EDTA.intact.gff3"),              emit: intact_gff3
        tuple val(meta), path("${meta}.renamed.ids.EDTA.pass.list"),    emit: pass_list
        tuple val(meta), path("${meta}.renamed.ids.EDTA.out"),          emit: out_file
        tuple val(meta), path("${meta}.EDTA.TElib.fa"),                 emit: te_list_fasta
        tuple val(meta), path("${meta}.renamed.ids.tsv"),               emit: renamed_ids_tsv

    script:
        """
        cat $pass_list > "${meta}.renamed.ids.EDTA.pass.list"
        cat $out_file > "${meta}.renamed.ids.EDTA.out"
        cat $te_lib_fa > "${meta}.EDTA.TElib.fa"
        cat $renamed_ids_tsv > "${meta}.renamed.ids.tsv"
        
        renamed_ids_head=\$(head -n 1 "$renamed_ids_tsv")
        
        if [[ \$renamed_ids_head == "IDs have acceptable length and character. No change required." ]]; then
            cat $te_anno_gff3 > "${meta}.EDTA.TEanno.gff3"
            cat $intact_gff3 > "${meta}.EDTA.intact.gff3"
        else
            reverse_edta_naming_f1b7bce.py "$renamed_ids_tsv" "$te_anno_gff3" "$intact_gff3" "$meta"
        fi
        """
}