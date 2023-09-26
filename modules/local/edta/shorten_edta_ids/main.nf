nextflow.enable.dsl=2

// https://github.com/Plant-Food-Research-Open/assembly_qc
// GPL-3.0: https://github.com/Plant-Food-Research-Open/assembly_qc/blob/main/LICENSE
process SHORTEN_EDTA_IDS {
    tag "${meta}"
    label "process_single"

    container "docker://gallvp/python3npkgs:v0.4"
    
    input:
        tuple val(meta), path(fasta_file)
    
    output:
        tuple val(meta), path("*.renamed.ids.fa"),  emit: renamed_ids_fasta
        tuple val(meta), path("*.renamed.ids.tsv"), emit: renamed_ids_tsv
    
    script:
        """
        fasta_file_bash_var="$fasta_file"
        output_prefix="\${fasta_file_bash_var%.*}"

        shorten_fasta_ids_c97537f.py "$fasta_file" "\$output_prefix"
        """
}