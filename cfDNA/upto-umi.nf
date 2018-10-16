#!/usr/bin/env nextflow
// Prerequisite: Build the genome index like so:
// /wittelab/software/hisat2-2.0.5/hisat2-build /wittelab/hg19/ucsc.hg19.fasta /wittelab/hg19/ucsc.hg19

import groovy.io.FileType

// Define the input channel
//Channel.fromFilePairs("${params.fastq_home}/*_R{1,2}*.fastq"){ file -> file.name.split('_')[0] }.println { id, files -> "Files with the sample ID $id are $files" }
read_pairs = Channel.fromFilePairs("${params.fastq_home}/${params.read_glob}"){ file -> file.name.split('_')[0] }
println("Using read data matching glob: ${params.fastq_home}/*_R{1,2}*.fastq.gz")


/************************
 * Alignment 
 ************************/
process Align {
  tag {"${sample_id}: ${reads[0].getName()} & ${reads[1].getName()}"}

  input:
  set sample_id, reads from read_pairs

  output:
  set sample_id, file("${sample_id}-sorted.bam") into aligned_files 

  """
  mkdir -p ${params.alignment_dir}
  mkdir -p ${params.alignment_dir}/logs
  export PERL5LIB=/wittelab/data1/software/lib/perl5
  export PU=`zcat < ${reads[0]} | head -n1 | awk 'BEGIN{FS=":"; OFS="."} /1/{print \$3,\$4,\$10}'`
  
  ${params.hisat2_home}/hisat2 \
    -q \
    -p 4 \
    -x ${params.hisat2_index} \
    -1 ${reads[0]} \
    -2 ${reads[1]} \
    --rg-id=${reads[0].getName()} \
    --rg PU:\$PU \
    --rg SM:${sample_id} \
    --rg PL:ILLUMINA \
    --rg LB:${sample_id} \
    2> ${params.alignment_dir}/logs/${reads[1].getName()}.log \
  | samtools view -Sb - | samtools sort -@ 4 -o ${sample_id}-sorted.bam
  """
}
// Group files by their sample ids 
joint_bams = aligned_files.groupTuple()


process MergeAligned {
  publishDir "${params.alignment_dir}", mode: "copy"
  tag {"${sample_id}"}

  input:
  set sample_id, bams from joint_bams

  output:
  set sample_id, file("${sample_id}-merged.bam"), file("${sample_id}-merged.bai") into merged_aligned

  """
  java -Xms4g -XX:ParallelGCThreads=4 -jar ${params.picard_home}/picard.jar MergeSamFiles \
    ASSUME_SORTED=true \
    CREATE_INDEX=true \
    INPUT=${bams.join(" INPUT=")} \
    MERGE_SEQUENCE_DICTIONARIES=false \
    OUTPUT=${sample_id}-merged.bam \
    SORT_ORDER=coordinate \
    USE_THREADING=true \
    VALIDATION_STRINGENCY=STRICT
  """
}


// Demultiplex reads by UMI and build consensus
process ProcessUMIs {
  publishDir "${params.alignment_dir}", mode: "copy"
  tag {"${sample_id}"}

  input:
  set sample_id, file(sample), file(bai) from merged_aligned

  output:
  set val("${sample_id - params.patient_regex}"), file("${sample_id}-deduped.bam") into deduped

  """
  mkdir -p ${params.alignment_dir}
  ${params.connor_home}/connor \
  -f ${params.consensus_freq_threshold} \
  -s ${params.min_family_size_threshold} \
  -d ${params.umt_distance_threshold} \
  ${sample} ${sample_id}-deduped.bam --force
  """
}
