# PAN-GENE
A NextFlow pipeline for pan-genome annotation.

## Pipeline Flowchart

```mermaid
flowchart LR
    TARGET_ASSEMBLIES((target_assemblies))
    TE_LIBRARIES((te_libraries))
    SAMPLESHEET((samplesheet))
    GUNZIP_TE[GUNZIP]
    SKIP_EDTA{Skip EDTA}
    pend((dev))
    
    TARGET_ASSEMBLIES --> GUNZIP
    GUNZIP --> FASTA_VALIDATE
    FASTA_VALIDATE --> FASTA_PERFORM_EDTA
    FASTA_VALIDATE --> SKIP_EDTA
    
    TE_LIBRARIES --> GUNZIP_TE
    GUNZIP_TE --> SKIP_EDTA
    SKIP_EDTA --> REPEATMASKER
    FASTA_PERFORM_EDTA --> REPEATMASKER
    REPEATMASKER --> BRAKER3

    SAMPLESHEET --> SAMPLESHEET_CHECK
    SAMPLESHEET_CHECK --> |Technical replicates|CAT_FASTQ
    CAT_FASTQ --> FASTQC
    SAMPLESHEET_CHECK --> FASTQC
    FASTQC --> FASTP
    FASTP --> pend

    subgraph Params
    TARGET_ASSEMBLIES
    TE_LIBRARIES
    SAMPLESHEET
    end

    subgraph Validate
    GUNZIP
    FASTA_VALIDATE
    end

    subgraph Repeatmask
    GUNZIP_TE
    FASTA_PERFORM_EDTA
    SKIP_EDTA
    REPEATMASKER
    end

    subgraph SamplePrep
    SAMPLESHEET_CHECK
    CAT_FASTQ
    FASTQC
    FASTP
    end

    style Params fill:#00FFFF21,stroke:#00FFFF21
    style Validate fill:#00FFFF21,stroke:#00FFFF21
    style Repeatmask fill:#00FFFF21,stroke:#00FFFF21
    style SamplePrep fill:#00FFFF21,stroke:#00FFFF21
```

## Plant&Food Users

Configure the pipeline by modifying `nextflow.config` and submit to SLURM for execution.

```bash
sbatch ./pan_gene_pfr.sh 
```


## Third-party Sources

Some software components of this pipeline have been adopted from following third-party sources:

1. nf-core [MIT](https://github.com/nf-core/modules/blob/master/LICENSE): https://github.com/nf-core/modules

> **The nf-core framework for community-curated bioinformatics pipelines.**
>
> Philip Ewels, Alexander Peltzer, Sven Fillinger, Harshil Patel, Johannes Alneberg, Andreas Wilm, Maxime Ulysse Garcia, Paolo Di Tommaso & Sven Nahnsen.
>
> _Nat Biotechnol._ 2020 Feb 13. doi: [10.1038/s41587-020-0439-x](https://dx.doi.org/10.1038/s41587-020-0439-x).

2. rewarewaannotation [MIT](https://github.com/kherronism/rewarewaannotation/blob/master/LICENSE): https://github.com/kherronism/rewarewaannotation
3. assembly_qc [GPL-3.0](https://github.com/Plant-Food-Research-Open/assembly_qc/blob/main/LICENSE): https://github.com/Plant-Food-Research-Open/assembly_qc