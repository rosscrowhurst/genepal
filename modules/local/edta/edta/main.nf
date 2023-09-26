nextflow.enable.dsl=2

// https://github.com/Plant-Food-Research-Open/assembly_qc
// GPL-3.0: https://github.com/Plant-Food-Research-Open/assembly_qc/blob/main/LICENSE
process EDTA {
    tag "$meta"
    label "process_high"
    label "process_week_long"
    
    container 'quay.io/biocontainers/edta:2.1.0--hdfd78af_1'
    containerOptions "-B $TMPDIR:$TMPDIR"

    input:
        tuple val(meta), path(fasta_file)
    
    output:
        tuple val(meta), path('*.EDTA.TEanno.gff3'),    emit: te_anno_gff3
        tuple val(meta), path('*.EDTA.intact.gff3'),    emit: intact_gff3
        tuple val(meta), path('*.EDTA.pass.list'),      emit: pass_list
        tuple val(meta), path('*.EDTA.out'),            emit: out_file
        tuple val(meta), path('*.EDTA.TElib.fa'),       emit: te_lib_fasta
    
    script:
        def args = task.ext.args ?: ''
        """
        EDTA.pl \\
        --genome $fasta_file \\
        --threads ${task.cpus} \\
        $args

        fasta_file_mod="${fasta_file}.mod"
        
        [[ -f "./\${fasta_file_mod}.EDTA.raw/LTR/\${fasta_file_mod}.pass.list" ]] \
        && echo "EDTA pass list detected" \
        || echo "EDTA PASS LIST IS EMPTY" > "./\${fasta_file_mod}.EDTA.raw/LTR/\${fasta_file_mod}.pass.list"
        
        ln -s "./\${fasta_file_mod}.EDTA.raw/LTR/\${fasta_file_mod}.pass.list" "\${fasta_file_mod}.EDTA.pass.list"
        
        ln -s "./\${fasta_file_mod}.EDTA.anno/\${fasta_file_mod}.out" "\${fasta_file_mod}.EDTA.out"
        """
}