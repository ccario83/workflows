#!/usr/bin/env nextflow
// Prerequisite: Build the genome index like so:
// /wittelab/software/hisat2-2.0.5/hisat2-build /wittelab/hg19/ucsc.hg19.fasta /wittelab/hg19/ucsc.hg19

import groovy.io.FileType

// Define the input channel of all read paired samples
//Channel.fromFilePairs("${params.fastq_home}/*_R{1,2}*.fastq"){ file -> file.name.split('_')[0] }.println { id, files -> "Files with the sample ID $id are $files" }
read_pairs = Channel.fromFilePairs("${params.fastq_home}/${params.read_glob}"){ file -> file.name.split('_')[0] }
//println("Using read data matching glob: ${params.fastq_home}/*_R{1,2}*.fastq.gz")


//************************
//* Alignment
//************************
process Align {
  tag {"${sample_id}: ${reads[0].getName()} & ${reads[1].getName()}"}
  
  input:
  set sample_id, reads from read_pairs

  output:
  set sample_id, file("${sample_id}.bam") into aligned_files 

  """  
  $SPEEDSEQ align \
    -o ${sample_id} \
    -t 32 \
    -R "@RG\tID:${reads[0].getName()}\tSM:${sample_id}\tLB:${sample_id}" \
    ${params.hg19_reference} \
    ${reads[0]} \
    ${reads[1]}
  """
}
// Group files by their sample ids 
joint_bams = aligned_files.groupTuple()

process MergeAligned {
  publishDir "${params.alignment_dir}", mode: "copy", pattern: "*.bam *.bai"
  tag {"${sample_id}"}

  input:
  set sample_id, bams from joint_bams

  output:
  set sample_id, file("${sample_id}_merged.bam"), file("${sample_id}_merged.bam.bai") into merged_aligned

  script:
  if(bams.size()==1)
    """
    mkdir -p ${params.alignment_dir}

    mv ${bams[0]} ${sample_id}_merged.bam 
    $SAMBAMBA index ${sample_id}_merged.bam
    """
  else
    """
    mkdir -p ${params.alignment_dir}

    $SAMBAMBA merge ${sample_id}_merged.bam ${bams.each{it[0]}.join(" ")}
    $SAMBAMBA index ${sample_id}_merged.bam
    """
}

normals = Channel.create()
tumors  = Channel.create()
paired  = Channel.create()

merged_aligned.choice( normals, tumors ) { a -> a =~ params.normal_regex ? 0 : 1 }
normals.cross(tumors){ a -> a[0].split("-")[0] }.into(paired)


//************************
//* Variant calling 
//************************
process CallVariants {
  publishDir "${params.variant_dir}", mode: "copy", pattern: "*.vcf.gz"
  tag {"${tumor[0]}"}

  input:
  set normal, tumor from paired

  output:
  file("${tumor[0]}.vcf.gz")

  """
  mkdir -p ${params.variant_dir}

  $SPEEDSEQ somatic \
    -o ${tumor[0]} \
    -F 0.05 \
    -q 1 \
    -t 32 \
    ${params.hg19_reference} \
    ${normal[1]} \
    ${tumor[1]} \
  """
}

