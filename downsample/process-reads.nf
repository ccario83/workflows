#!/usr/bin/env nextflow
import groovy.io.FileType


read_pairs = Channel.fromFilePairs("${params.fastq_home}/*_R{1,2}.fastq.gz"){file -> file.name.split('-')[0]}
println("Using read data matching glob: ${params.fastq_home}/*_R{1,2}.fastq.gz")

/************************
 * Downsampling
 ************************/
process Downsampling {
  tag {"${sample_id} => ${reads.collect {it.getName()}}"}
  storeDir params.output_dir

  input:
  set sample_id, reads from read_pairs

  output:
  set file("${sample_id}_R1.fastq.gz"), file("${sample_id}_R2.fastq.gz")
 
  """
  desired=\$(cat ${params.sample_reads} | grep ${sample_id} | sed 's/^.*\\t\\(.*\\)\$/\\1/')
  seqtk sample -s101 ${reads[0]} \$desired | gzip > ${sample_id}_R1.fastq.gz
  seqtk sample -s101 ${reads[1]} \$desired | gzip > ${sample_id}_R2.fastq.gz
  """
}
