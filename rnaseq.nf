#!/usr/bin/env nextflow
VERSION="1.0"
nextflow.enable.dsl=2

// Include your modules and parameterize them
include { FASTQC }                from './modules/FASTQC' params(params)
include { TRIM_GALORE }           from './modules/TRIM_GALORE' params(params)
include { HISAT_INDEX }           from './modules/HISAT' params(params)
include { HISAT_MAPPING }         from './modules/HISAT' params(params)
include { FEATURE_COUNTS }        from './modules/QUANT' params(params)
include { FEATURE_COUNTS_MATRIX } from './modules/QUANT' params(params)
include { ESET }                  from './modules/QUANT' params(params)
include { MULTIQC }               from './modules/MULTIQC' params(params)

//workflow {
//  HISAT_INDEX( params.fasta, params.gtf )
//}

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

// paired-end reads
workflow {
  load_reads("${params.wd}/raw_data/rnaseq/reads/*_{1,2}.fq.gz", params.paired)
  TRIM_GALORE( reads )
  HISAT_MAPPING( TRIM_GALORE.out[0] )
  FEATURE_COUNTS( HISAT_MAPPING.out[0], params.gtf)
  FEATURE_COUNTS_MATRIX( FEATURE_COUNTS.out[0].collect() )
  ESET( FEATURE_COUNTS_MATRIX.out[0] )
  MULTIQC( ESET.out[0] )
}

/*
// single-end reads
workflow {
  load_reads("${params.wd}/data/rnaseq/reads/*_1.fq.gz", params.paired)
  TRIM_GALORE( reads )
  HISAT_MAPPING( TRIM_GALORE.out[0] )
  FEATURE_COUNTS( HISAT_MAPPING.out[0], params.gtf)
  FEATURE_COUNTS_MATRIX( FEATURE_COUNTS.out[0].collect() )
  ESET( FEATURE_COUNTS_MATRIX.out[0] )
  MULTIQC( ESET.out[0] )
}
*/

/*
// start from bam files
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
*/
