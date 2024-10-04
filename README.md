[![GitHub Actions CI Status](https://github.com/plant-food-research-open/genepal/actions/workflows/ci.yml/badge.svg)](https://github.com/plant-food-research-open/genepal/actions/workflows/ci.yml)
[![GitHub Actions Linting Status](https://github.com/plant-food-research-open/genepal/actions/workflows/linting.yml/badge.svg)](https://github.com/plant-food-research-open/genepal/actions/workflows/linting.yml)[![Cite with Zenodo](http://img.shields.io/badge/DOI-10.5281/zenodo.XXXXXXX-1073c8?labelColor=000000)](https://doi.org/10.5281/zenodo.XXXXXXX)
[![nf-test](https://img.shields.io/badge/unit_tests-nf--test-337ab7.svg)](https://www.nf-test.com)

[![Nextflow](https://img.shields.io/badge/nextflow%20DSL2-%E2%89%A523.04.0-23aa62.svg)](https://www.nextflow.io/)
[![run with conda ❌](http://img.shields.io/badge/run%20with-conda%20❌-3EB049?labelColor=000000&logo=anaconda)](https://docs.conda.io/en/latest/)
[![run with docker](https://img.shields.io/badge/run%20with-docker-0db7ed?labelColor=000000&logo=docker)](https://www.docker.com/)
[![run with singularity](https://img.shields.io/badge/run%20with-singularity-1d355c.svg?labelColor=000000)](https://sylabs.io/docs/)
[![Launch on Seqera Platform](https://img.shields.io/badge/Launch%20%F0%9F%9A%80-Seqera%20Platform-%234256e7)](https://cloud.seqera.io/launch?pipeline=https://github.com/plant-food-research-open/genepal)

## Introduction

**plant-food-research-open/genepal** is a bioinformatics pipeline for single genome, multiple genomes and pan-genome annotation. An overview is shown in the [Pipeline Flowchart](#pipeline-flowchart) and the references for the tools are listed in [CITATIONS.md](./CITATIONS.md).

## Pipeline Flowchart

<p align="center"><img src="docs/img/genepal.png"></p>

- [FASTA VALIDATOR](https://github.com/linsalrob/fasta_validator): Validate genome fasta
- [REPEATMODELER](https://github.com/Dfam-consortium/RepeatModeler) or [EDTA](https://github.com/oushujun/EDTA): Create TE library
- [REPEATMASKER](https://github.com/rmhubley/RepeatMasker): Soft mask the genome fasta
- [FASTQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc), [FASTP](https://github.com/OpenGene/fastp), [SORTMERNA](https://github.com/sortmerna/sortmerna): QC, trim and filter RNASeq evidence
- [STAR](https://github.com/alexdobin/STAR): RNASeq alignment
- [BRAKER](https://github.com/Gaius-Augustus/BRAKER): Annotate the genome fasta
- [LIFTOFF](https://github.com/agshumate/Liftoff): Liftoff annotations from reference genome fasta/gff
- [TSEBRA](https://github.com/Gaius-Augustus/TSEBRA), [AGAT](https://github.com/NBISweden/AGAT): Merge BRAKER and Liftoff annotations
- [EGGNOG-MAPPER](https://github.com/eggnogdb/eggnog-mapper): Add functional annotation to gff
- [ORTHOFINDER](https://github.com/davidemms/OrthoFinder): Perform phylogenetic orthology inference across input genomes
- [GENOMETOOLS](https://github.com/genometools/genometools), [GFFREAD](https://github.com/gpertea/gffread): Final GFF format validation and extraction of protein sequences
- [BUSCO](https://gitlab.com/ezlab/busco): Completeness statistics for genome and annotation through proteins

## Usage

Refer to [usage](./docs/usage.md), [parameters](./docs/parameters.md) and [output](./docs/output.md) documents for details.

> [!NOTE]
> If you are new to Nextflow and nf-core, please refer to [this page](https://nf-co.re/docs/usage/installation) on how to set-up Nextflow. Make sure to [test your setup](https://nf-co.re/docs/usage/introduction#how-to-run-a-pipeline) with `-profile test` before running the workflow on actual data.

First, prepare an assemblysheet with your input genomes that looks as follows:

`assemblysheet.csv`:

```csv
tag         ,fasta              ,is_masked
a_thaliana  ,/path/to/genome.fa ,yes
```

Each row represents an input genome and the fields are:

- `tag:` A unique tag which represents the genome throughout the pipeline
- `fasta:` fasta file for the genome
- `is_masked`: yes or no to denote whether the fasta file is already masked or not

At minimum, a file with proteins as evidence is also required. Now, you can run the pipeline using:

```bash
nextflow run plant-food-research-open/genepal \
  -profile <docker/singularity/.../institute> \
  --input assemblysheet.csv \
  --protein_evidence proteins.faa \
  --outdir <OUTDIR>
```

> [!WARNING]
> Please provide pipeline parameters via the CLI or Nextflow `-params-file` option. Custom config files including those provided by the `-c` Nextflow option can be used to provide any configuration _**except for parameters**_;
> see [docs](https://nf-co.re/usage/configuration#custom-configuration-files).

### Plant&Food Users

Download the pipeline to your `/workspace/$USER` folder. Change the parameters defined in the [pfr/params.json](./pfr/params.json) file. Submit the pipeline to SLURM for execution.

```bash
sbatch ./pfr_genepal
```

## Credits

plant-food-research-open/genepal workflows were originally scripted by Jason Shiller ([@jasonshiller](https://github.com/jasonshiller)). Usman Rashid ([@gallvp](https://github.com/gallvp)) wrote the Nextflow pipeline.

We thank the following people for their extensive assistance in the development of this pipeline:

- Cecilia Deng [@CeciliaDeng](https://github.com/CeciliaDeng)
- Charles David [@charlesdavid](https://github.com/charlesdavid)
- Chen Wu [@christinawu2008](https://github.com/christinawu2008)
- Leonardo Salgado [@leorippel](https://github.com/leorippel)
- Ross Crowhurst [@rosscrowhurst](https://github.com/rosscrowhurst)
- Susan Thomson [@cflsjt](https://github.com/cflsjt)
- Ting-Hsuan Chen [@ting-hsuan-chen](https://github.com/ting-hsuan-chen)

The pipeline uses nf-core modules contributed by following authors:

<a href="https://github.com/gallvp"><img src="https://github.com/gallvp.png" width="50" height="50"></a>
<a href="https://github.com/drpatelh"><img src="https://github.com/drpatelh.png" width="50" height="50"></a>
<a href="https://github.com/kevinmenden"><img src="https://github.com/kevinmenden.png" width="50" height="50"></a>
<a href="https://github.com/adamrtalbot"><img src="https://github.com/adamrtalbot.png" width="50" height="50"></a>
<a href="https://github.com/toniher"><img src="https://github.com/toniher.png" width="50" height="50"></a>
<a href="https://github.com/joseespinosa"><img src="https://github.com/joseespinosa.png" width="50" height="50"></a>
<a href="https://github.com/grst"><img src="https://github.com/grst.png" width="50" height="50"></a>
<a href="https://github.com/edmundmiller"><img src="https://github.com/edmundmiller.png" width="50" height="50"></a>
<a href="https://github.com/maxulysse"><img src="https://github.com/maxulysse.png" width="50" height="50"></a>
<a href="https://github.com/kherronism"><img src="https://github.com/kherronism.png" width="50" height="50"></a>
<a href="https://github.com/vagkaratzas"><img src="https://github.com/vagkaratzas.png" width="50" height="50"></a>
<a href="https://github.com/robsyme"><img src="https://github.com/robsyme.png" width="50" height="50"></a>
<a href="https://github.com/priyanka-surana"><img src="https://github.com/priyanka-surana.png" width="50" height="50"></a>
<a href="https://github.com/praveenraj2018"><img src="https://github.com/praveenraj2018.png" width="50" height="50"></a>
<a href="https://github.com/muffato"><img src="https://github.com/muffato.png" width="50" height="50"></a>
<a href="https://github.com/matthdsm"><img src="https://github.com/matthdsm.png" width="50" height="50"></a>
<a href="https://github.com/mashehu"><img src="https://github.com/mashehu.png" width="50" height="50"></a>
<a href="https://github.com/mahesh-panchal"><img src="https://github.com/mahesh-panchal.png" width="50" height="50"></a>
<a href="https://github.com/jvhagey"><img src="https://github.com/jvhagey.png" width="50" height="50"></a>
<a href="https://github.com/jfy133"><img src="https://github.com/jfy133.png" width="50" height="50"></a>
<a href="https://github.com/jemten"><img src="https://github.com/jemten.png" width="50" height="50"></a>
<a href="https://github.com/friederikehanssen"><img src="https://github.com/friederikehanssen.png" width="50" height="50"></a>
<a href="https://github.com/felixkrueger"><img src="https://github.com/felixkrueger.png" width="50" height="50"></a>
<a href="https://github.com/ewels"><img src="https://github.com/ewels.png" width="50" height="50"></a>
<a href="https://github.com/erikrikarddaniel"><img src="https://github.com/erikrikarddaniel.png" width="50" height="50"></a>
<a href="https://github.com/charles-plessy"><img src="https://github.com/charles-plessy.png" width="50" height="50"></a>

## Contributions and Support

If you would like to contribute to this pipeline, please see the [contributing guidelines](.github/CONTRIBUTING.md).

## Citations

<!-- TODO nf-core: Add citation for pipeline after first release. Uncomment lines below and update Zenodo doi and badge at the top of this file. -->
<!-- If you use plant-food-research-open/genepal for your analysis, please cite it using the following doi: [10.5281/zenodo.XXXXXX](https://doi.org/10.5281/zenodo.XXXXXX) -->

An extensive list of references for the tools used by the pipeline can be found in the [`CITATIONS.md`](CITATIONS.md) file.

This pipeline uses code and infrastructure developed and maintained by the [nf-core](https://nf-co.re) community, reused here under the [MIT license](https://github.com/nf-core/tools/blob/master/LICENSE).

> **The nf-core framework for community-curated bioinformatics pipelines.**
>
> Philip Ewels, Alexander Peltzer, Sven Fillinger, Harshil Patel, Johannes Alneberg, Andreas Wilm, Maxime Ulysse Garcia, Paolo Di Tommaso & Sven Nahnsen.
>
> _Nat Biotechnol._ 2020 Feb 13. doi: [10.1038/s41587-020-0439-x](https://dx.doi.org/10.1038/s41587-020-0439-x).
