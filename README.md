# Pipeliner DSL2

<i>[Pipeliner](https://github.com/montilab/pipeliner) upgraded for Nextflow DSL2 modules</i>   

[![Built With](https://img.shields.io/badge/nextflow-DSL2-23aa62.svg?labelColor=000000)](https://www.nextflow.io/)
![Compatibility](https://img.shields.io/badge/Compatibility-Linux%20%2F%20OSX-orange.svg)
![Dependencies](https://img.shields.io/badge/dependencies-up%20to%20date-brightgreen.svg)
[![GitHub Issues](https://img.shields.io/github/issues/montilab/pipeliner-2.svg)](https://github.com/montilab/pipeliner-2/issues)

# Overview

We developed the _yQTL Pipeline_ – with ‘y’ indicating the dependent quantitative variable being modeled – to facilitate and automate large-scale QTL analysis. Prior to genome-wide association test, the pipeline supports the calculation or the direct input of pre-defined genome-wide principal components and genetic relationship matrix when applicable. User-specified covariates can also be provided. Depending on whether familial relatedness exists among the subjects or not, genome-wide association tests will be performed using either a linear mixed-effect model or a linear model, respectively. Through the adoption of the workflow management tool Nextflow, the pipeline parallelizes the analysis steps to optimize run-time and ensure reproducibility of the results. A user-friendly R Shiny App is also provided for the visualization of the results, including Manhattan plots of user-selected phenotype traits, and trait-QTL connection networks based on user-specified p-value thresholds.

## Clone Repository

```bash
$ git clone https://github.com/montilab/pipeliner-2
```

## Installing Nextflow

Workflows are built using [Nextflow](https://www.nextflow.io/). Nextflow can be used on any POSIX compatible system (Linux, OS X, etc) and requires BASH and Java 8 (or higher) to be installed. Download the latest version of Nextflow compatible with DSL2:

1. Make sure 8 or later is installed on your computer by using the command.

   ```bash
   java -version
   ```

2. Download the nextflow executable, this repository uses nextflow version 20.10.0.5430.

   ```bash
   curl -s https://get.nextflow.io | bash
   ```

3. Add nextflow to your $PATH in your .bash_profile (or system equivalent file).

   ```bash
   export PATH=$PATH:/path/to/nextflow
   ```

## Nextflow DSL2

To make use of this repository you'll need a basic understanding of how [processes and channels](https://www.nextflow.io/docs/latest/basic.html) work in Nextflow as well as familiarity with the new syntax for defining modular workflows [Nextflow DSL2](https://www.nextflow.io/docs/latest/dsl2.html). Also note that many of the modules will contain Groovy code snippets.

## Modular Workflows

A Nextflow workflow is a series of dependent processes that execute pre-defined scripts. Before DSL2, we wrote pipelines that shared script templates to reuse software tools common across pipelines. Recently, it is now possible to write processes as modules (similar to functions). This allows us to make modular pipelines and great reduces the friction of modifying or creating new pipelines from existing components.

A couple important notes about the beginning of your workflow script...

1. You must have `nextflow.enable.dsl=2` in the beginning to enable DSL2.
2. Your modules will be parameterized with the `params` variable. So check the modules that you're using and see what parameters they are expecting to be defined.
3. The `params` variable must come before the module imports, because they're parameterized as they're imported.
4. Here I am importing them all, but you only need to import the modules you're making use of.

```groovy
#!/usr/bin/env nextflow
VERSION="1.0"
nextflow.enable.dsl=2

params.paired = false
params.wd = "path/to/pipeliner-2"
params.outdir = "${params.wd}/results"
params.fasta = "${params.wd}/data/genomes/genome_reference.fa"
params.gtf = "${params.wd}/data/genomes/genome_annotation.gtf"
params.index = "${params.wd}/results/hisat/index/part"

include { FASTQC }                from './modules/FASTQC' params(params)
include { TRIM_GALORE }           from './modules/TRIM_GALORE' params(params)
include { HISAT_INDEX }           from './modules/HISAT' params(params)
include { HISAT_MAPPING }         from './modules/HISAT' params(params)
include { FEATURE_COUNTS }        from './modules/QUANT' params(params)
include { FEATURE_COUNTS_MATRIX } from './modules/QUANT' params(params)
include { ESET }                  from './modules/QUANT' params(params)
include { MULTIQC }               from './modules/MULTIQC' params(params)
```

Here are some examples of how you might use these modules.

### Building an Index

You'll probably already have an index built but in case you need to build one, you can now do it independent of your workflow.

```groovy
workflow {
  HISAT_INDEX( params.fasta, params.gtf )
}
```

### Processing Paired RNA-seq Reads

We defined a helper function for reading in single or paired end reads. All you have to do is specify the files you want to read in and make sure `params.paired = true` and your reads will be properly formatted for the downstream modules.

First we pass the reads to `TRIM_GALORE` which trims off the adapters and performs some quality control checks. The first output are the trimmed reads so we can pass `TRIM_GALORE.out[0]` to `HISAT_MAPPING`. Check your modules to see what the inputs/outputs are so you can modify this sequence if necessary.

Note that we pass `FEATURE_COUNTS.out[0].collect()` to `FEATURE_COUNTS_MATRIX`. The `.collect()` function will flatten the counts per sample so that they are all passed at once and processed together rather than processed in independently in parallel. 

```groovy
def load_reads(path, paired) {
  if (paired) {
    Channel
      .fromFilePairs( path )
      .set {reads}
  } else {
    Channel
      .fromPath( path )
      .map { [it.getName().split("\\_1|\\_2", 2)[0], [it]] }
      .set {reads}
  }
}

workflow {
  load_reads("${params.wd}/data/rnaseq/reads/*_{1,2}.fq.gz", params.paired)
  TRIM_GALORE( reads )
  HISAT_MAPPING( TRIM_GALORE.out[0] )
  FEATURE_COUNTS( HISAT_MAPPING.out[0], params.gtf)
  FEATURE_COUNTS_MATRIX( FEATURE_COUNTS.out[0].collect() )
  ESET( FEATURE_COUNTS_MATRIX.out[0] )
  MULTIQC( ESET.out[0] )
}
```

If the reads were single end, we only have to change the wildcard and ensure `params.paired = false`.

```groovy
load_reads("${params.wd}/data/rnaseq/reads/*_1.fq.gz", params.paired)
```

### Process Just the Bam Files

What if we already mapped the reads and have the bams? Rather than writing this conditional logic into a pipeline explicity, with a modular architecture, we can just modify the modules we're using.

```groovy
def load_bams(path) {
  Channel
      .fromPath( path )
      .map { [it.getName().split("\\.bam", 0)[0], [it]] }
      .set {bams}
}

workflow {
  load_bams("${params.wd}/data/rnaseq/bams/*.bam")
  FEATURE_COUNTS( bams, params.gtf)
  FEATURE_COUNTS_MATRIX( FEATURE_COUNTS.out[0].collect() )
  ESET( FEATURE_COUNTS_MATRIX.out[0] )
}
```
