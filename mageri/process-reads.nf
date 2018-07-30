#!/usr/bin/env nextflow
read_pairs = Channel.fromFilePairs("${params.fastq_home}/${params.read_glob}"){ file -> file.name.split('_')[0] }

//lane_pairs = Channel.fromFilePairs("${params.fastq_home}/${params.lane_glob}", size: -1){ file -> file.name.split('_')[0] + "_" + file.name.split('_')[3]}
//println("Using read data matching glob: ${params.fastq_home}/${params.lane_glob}")
// Merging is not needed since combined version of all files now exist through the combined workflow.
// To include this process, channels must be modified.
/************************
 * Merging
 ************************
process Combine {
  tag {"Merging ${reads.collect {it.getName()}}"}
  clusterOptions "-l vmem=16gb,mem=16gb"

  input:
  set sample_id, reads from lane_pairs

  output:
  set sample_id, file("${sample_id}.fastq.gz") into joint_reads

  script:
  """
  #cat ${reads.join(" ")} > ${sample_id}.fastq.gz
  touch ${sample_id}.fastq.gz
  """
}
// Group files by their sample ids
grouped_reads = joint_reads.groupTuple()
*/

//************************
//* Mageri UMI variant caller
//************************
process Mageri {
  tag {"${sample_id}: ${reads[0].getName()} & ${reads[1].getName()}"}
  clusterOptions "-l vmem=200gb,mem=128gb,nodes=1:ppn=12"

  input:
  set sample_id, reads from read_pairs

  """
  java -Xms128G -Xmx200G -XX:+UseSerialGC -XX:CICompilerCount=4 -jar mageri.jar \
  -M3 ${params.barcode_pattern} \
  --import-preset ${params.preset} \
  --references ${params.panel_reference} \
  --bed ${params.panel_loci} \
  --project-name ${params.project_name} \
  --sample-name ${sample_id} \
  -R1 ${reads[0]} \
  -R2 ${reads[1]} \
  -O ${params.variant_dir}/
  """
}
