nextflow.enable.dsl=2

// https://github.com/Plant-Food-Research-Open/assembly_qc
// GPL-3.0: https://github.com/Plant-Food-Research-Open/assembly_qc/blob/main/LICENSE
process SHORTEN_EDTA_IDS {
    tag "$meta.id"
    label "process_single"

    container "docker://gallvp/python3npkgs:v0.4"
    
    input:
        tuple val(meta), path(fasta_file)
    
    output:
        tuple val(meta), path("*.renamed.ids.fa"),  emit: renamed_ids_fasta
        tuple val(meta), path("*.renamed.ids.tsv"), emit: renamed_ids_tsv
        path "versions.yml",                        emit: versions
    
    script:
        def VERSION = "c97537f" // WARN: Version information not provided by tool on CLI. Please update version string below when bumping container versions.
        """
        FILE="$fasta_file"
        output_prefix="\${FILE%.*}"

        shorten_fasta_ids_c97537f.py "$fasta_file" "\$output_prefix"

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            shorten_fasta_ids: $VERSION
        END_VERSIONS
        """
    
    stub:
        def VERSION = "c97537f" // WARN: Version information not provided by tool on CLI. Please update version string below when bumping container versions.
        """
        FILE="$fasta_file"
        output_prefix="\${FILE%.*}"

        touch "\${output_prefix}.renamed.ids.fa"
        touch "\${output_prefix}.renamed.ids.tsv"

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            shorten_fasta_ids: $VERSION
        END_VERSIONS
        """
}