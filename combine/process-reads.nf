#!/usr/bin/env nextflow
import groovy.io.FileType


lane_pairs = Channel.fromFilePairs("${params.fastq_home}/${params.read_glob}", size: -1){ file -> file.name.split('_')[0] + "_" + file.name.split('_')[3]}
println("Using read data matching glob: ${params.fastq_home}/${params.read_glob}")

/************************
 * Alignment 
 ************************/
process Combine {
  tag {"Merging ${reads.collect {it.getName()}}"}
  storeDir params.output_dir

  input:
  set sample_id, reads from lane_pairs

  output:
  file("${sample_id}.fastq.gz")

  script: 
  """
  cat ${reads.join(" ")} > ${sample_id}.fastq.gz
  """
}
