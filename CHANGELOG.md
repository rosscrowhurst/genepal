# PlantandFoodResearch/pangene: Changelog

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## 0.4 - [21-Jun-2024]

### `Added`

1. Added `orthofinder_annotations` param
2. Added `FASTA_GFF_ORTHOFINDER` sub-workflow
3. Updated modules: `AGAT_CONVERTSPGFF2GTF`, `CAT_FASTQ`, `CUSTOM_DUMPSOFTWAREVERSIONS`, `EGGNOGMAPPER`, `FASTP`, `GFFREAD`, `SAMTOOLS_CAT`, `CUSTOM_RESTOREGFFIDS`, `CUSTOM_SHORTENFASTAIDS` and `EDTA_EDTA`
4. Updated sub-workflows: `FASTQ_FASTQC_UMITOOLS_FASTP` and `FASTA_EDTA_LAI`
5. Added evaluation by BUSCO [#41](https://github.com/PlantandFoodResearch/pangene/issues/41)
6. Included common tax ids for eggnog mapper [#27](https://github.com/PlantandFoodResearch/pangene/issues/27)
7. Implemented hierarchical naming scheme: geneI.tJ, geneI.tJ.exonK, geneI.tJ.cdsK [#19](https://github.com/PlantandFoodResearch/pangene/issues/19), [#34](https://github.com/PlantandFoodResearch/pangene/issues/34)
8. Now sorting list of bam and list of fastq before cat to avoid resume cache misses
9. Added a local modification to `FASTQ_FASTQC_UMITOOLS_FASTP` for stub and PR to nf-core/modules <https://github.com/nf-core/modules/pull/5858>
10. Allowed BAM files for RNA evidence [#3](https://github.com/PlantandFoodResearch/pangene/issues/3)

### `Fixed`

1. Fixed BRAKER spellings [#36](https://github.com/PlantandFoodResearch/pangene/issues/36)
2. Fixed liftoff failure when lifting off from a single reference [#40](https://github.com/PlantandFoodResearch/pangene/issues/40)
3. Added versions from GFF_STORE sub-workflows [#33](https://github.com/PlantandFoodResearch/pangene/issues/33)

### `Dependencies`

1. NextFlow!>=23.04.4
2. nf-validation=1.1.3

### `Deprecated`

1. Renamed `external_protein_fastas` param to `protein_evidence`
2. Renamed `fastq` param to `rna_evidence`

## 0.3.3 - [18-Jun-2024]

### `Added`

1. Added a stub test to evaluate the case where an assembly is soft masked but has no annotations

### `Fixed`

1. Fixed a bug where `is_masked` was ignored by the pipeline
2. Fixed a bug in param validation which allowed specification of `braker_hints` without `braker_gff3`

### `Dependencies`

1. NextFlow!>=23.04.4
2. nf-validation=1.1.3

### `Deprecated`

## 0.3.2 - [13-May-2024]

### `Added`

### `Fixed`

1. Increased time limit for REPEATMODELER_REPEATMODELER to 5 days
2. Now removing comments from fasta file before feeding it to BRAKER added tests for the perl one liner
3. Fixed CHANGELOG version check failure in `version_check.sh`
4. Increased the SLURM job time limit to 14 days

### `Dependencies`

1. NextFlow!>=23.04.4
2. nf-validation=1.1.3

### `Deprecated`

## 0.3.1 - [10-May-2024]

### `Added`

### `Fixed`

1. Increased time limit for REPEATMODELER_REPEATMODELER to 3 days, REPEATMASKER to 2 days, EDTA_EDTA to 7 days, BRAKER3 to 7 days and EGGNOGMAPPER to 1 day

### `Dependencies`

1. NextFlow!>=23.04.4
2. nf-validation=1.1.3

### `Deprecated`

## 0.3.0 - [30-April-2024]

### `Added`

1. Added changelog and semantic versioning
2. Changed license to MIT
3. Updated `.editorconfig`
4. Moved .literature to test/ branch
5. Renamed `pangene_local` to `local_pangene`
6. Renamed `pangene_pfr` to `pfr_pangene`
7. Added versioning checking
8. Updated github workflow to use pre-commit instead of prettier and editorconfig check
9. Added central singularity cache dir for pfr config
10. Added `SORTMERNA_INDEX` before `SORTMERNA`
11. Fixed sample contamination bug introduced by `file.simpleName`
12. Now using empty files for stub testing in CI
13. Now BRAKER can be skipped by including BRAKER outputs from previous runs in the `target_assemblies` param
14. Added `gffcompare` to merge liftoff annotations
15. Renamed `samplesheet` param to `fastq`
16. Now using assemblysheet in combination with nf-validation for assembly input
17. Added nextflow_schema.json
18. Now using nf-validation to validate fastqsheet provided by params.fastq
19. Moved `manifest.config` and `reporting_defaults.config` content to `nextflow.config`
20. Now using a txt file for `params.external_protein_fastas`
21. Now using nf-validation for `params.liftoff_annotations`
22. Now using nf-validation for all the parameters
23. Added `PURGE_BRAKER_MODELS` sub-workflow
24. Added `GFF_EGGNOGMAPPER` sub-workflow
25. Now using a custom version of `GFFREAD` which supports `meta` and `fasta`
26. Now using TSEBRA to purge models which do not have full intron support from BRAKER hints
27. Added params `eggnogmapper_evalue` and `eggnogmapper_pident`
28. Added `PURGE_NOHIT_BRAKER_MODELS` sub-workflow
29. Now merging BRAKER and liftoff models before running eggnogmapper
30. Added `GFF_MERGE_CLEANUP` sub-workflow
31. Now using `description` field to store notes and textual annotations in the gff files
32. Now using `mRNA` in place of `transcript` in gff files
33. Now `eggnogmapper_purge_nohits` is set to `false` by default
34. Added `GFF_STORE` sub workflow
35. `external_protein_fastas` and `eggnogmapper_db_dir` are not mandatory parameters
36. Added contributors
37. Add a document for the pipeline parameters
38. Updated `pfr_pangene` and `pfr/profile.config`
39. Now using local tests/stub files for GitHub CI
40. Now removing iso-forms left by TSEBRA using `AGAT_SPFILTERFEATUREFROMKILLLIST`
41. Added `pyproject.toml`
42. Now using PFAMs from eggnog if description is '-'

### `Fixed`

1. Removed liftoff models with `valid_ORF=False`
2. Updated license text to include 'Copyright (c) 2024 The New Zealand Institute for Plant and Food Research Limited'

### `Dependencies`

1. NextFlow!>=23.04.4
2. nf-validation=1.1.3

### `Deprecated`
