#!/usr/bin/env nextflow
read_pairs = Channel.fromFilePairs("${params.fastq_home}/${params.read_glob}"){ file -> file.name.split('_')[0] }

//************************
//* Alignment
//************************
process SpeedAlign {
  tag {"${sample_id}: ${reads[0].getName()} & ${reads[1].getName()}"}
  clusterOptions "-l vmem=64gb,mem=64gb,nodes=1:ppn=12"
  publishDir "${params.alignment_dir}", mode: "copy"
  
  input:
  set sample_id, reads from read_pairs

  output:
  set sample_id, file("${sample_id}.bam"), file("${sample_id}.bam.bai") into aligned 

  """  
  $SPEEDSEQ align \
    -o ${sample_id} \
    -t 10 \
    -R "@RG\tID:${reads[0].getName()}\tSM:${sample_id}\tLB:${sample_id}" \
    ${params.hg19_reference} \
    ${reads[0]} \
    ${reads[1]}
  """
}


//************************
//* Collapse UMI families
//************************
process ProcessUMIs {
  tag {"${sample_id}"}
  clusterOptions "-l vmem=64gb,mem=64gb"
  publishDir "${params.collapsed_dir}", mode: "copy"

  input:
  set sample_id, file(sample), file(sample_index) from aligned

  output:
  set sample_id, file("${sample_id}-collapsed.bam"), file("${sample_id}-collapsed.bam.bai")  into collapsed

  """
  mkdir -p ${params.collapsed_dir}
  connor \
    -f ${params.consensus_freq_threshold}\
    -s ${params.min_family_size_threshold}\
    -d {params.umt_distance_threshold}\
    --force \
    ${sample} \
    ${sample_id}-collapsed.bam 
  """
}


//************************
//* Variant calling 
//************************
process Call {
  tag {"${sample_id}"}
  clusterOptions "-l vmem=256gb,mem=256gb"
  publishDir "${params.variant_dir}", mode: "copy", pattern: "*.vcf.gz"
  

  input:
  set sample_id, file(sample), file(sample_index) from collapsed

  output:
  set sample_id, file("${sample_id}.vcf")

  """
  mkdir -p ${params.variant_dir}
  
  freebayes-1.0\
    -f ${params.hg19_reference}\
    -F ${params.min_alternate_fraction}\
    -C ${params.min_alternate_count}\
    ${sample} > ${sample_id}.vcf
  """
}


/*
// NOTE: Does not support --min_alternate_fraction and --min_alternate_count arguments
//************************
//* Variant calling 
//************************
process SpeedCall {
  tag {"${tumor[0]}"}
  clusterOptions "-l vmem=64gb,mem=64gb,nodes=1:ppn=32"
  publishDir "${params.variant_dir}", mode: "copy", pattern: "*.vcf.gz"

  input:
  set sample_id, file(sample), file(sample_index) from collapsed

  output:
  file("${sample_id}.vcf.gz")

  """
  mkdir -p ${params.variant_dir}
  export LD_LIBRARY_PATH="/opt/gcc/gcc-4.8.1/lib64/"
  
  $SPEEDSEQ var \
    -o ${sample_id} \
    -q 1 \
    -t 30 \
    ${params.hg19_reference} \
    ${sample}
  """
}
*/
