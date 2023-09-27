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
        tuple val(meta), path("$validFasta"),   emit: valid_fasta
        path "versions.yml",                    emit: versions

    script:
        validFasta = (fasta_file.toString() - ~/\.\w+$/) + ".validated.fasta"
        def VERSION = "a6a2ec1_ps" // WARN: Version information not provided by tool on CLI. Please update version string below when bumping container versions.
        """
        fasta_validate -v $fasta_file >/dev/null

        # If invalid, the above command will fail and
        # the NXF error startegy will kick in.
        
        cat $fasta_file > $validFasta

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            fasta_validate: $VERSION
        END_VERSIONS
        """
    
    stub:
        validFasta = (fasta_file.toString() - ~/\.\w+$/) + ".validated.fasta"
        def VERSION = "a6a2ec1_ps" // WARN: Version information not provided by tool on CLI. Please update version string below when bumping container versions.
        """
        touch $validFasta

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            fasta_validate: $VERSION
        END_VERSIONS
        """
}