nextflow.enable.dsl=2

// https://github.com/Plant-Food-Research-Open/assembly_qc
// GPL-3.0: https://github.com/Plant-Food-Research-Open/assembly_qc/blob/main/LICENSE
process FASTA_VALIDATE {
    tag "$meta.id"
    label "process_single"

    container "docker://gallvp/fasta_validator:a6a2ec1_ps"

    input:
        tuple val(meta), path(fasta_file)
    
    output:
        val(meta)

    script:
        """
        fasta_validate -v $fasta_file >/dev/null

        # If invalid, the above command will fail and
        # the NXF error startegy will kick in.
        
        echo -n "VALIDATE_FASTA:${meta.id}:VALID"
        """
}