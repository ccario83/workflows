#!/usr/bin/env nextflow
// Prerequisite: Build the genome index like so:
// /wittelab/software/hisat2-2.0.5/hisat2-build /wittelab/hg19/ucsc.hg19.fasta /wittelab/hg19/ucsc.hg19

import groovy.io.FileType

// Get a list of fastq files
def fastqs = []
def dir = new File("${params.fastq_home}").eachFileRecurse (FileType.FILES) { file -> fastqs << file.name.split('_')[0] }
def normal_lookup = match(fastqs.unique().collect())


// Define the input channel of all read paired samples
//Channel.fromFilePairs("${params.fastq_home}/*_R{1,2}*.fastq"){ file -> file.name.split('_')[0] }.println { id, files -> "Files with the sample ID $id are $files" }
read_pairs = Channel.fromFilePairs("${params.fastq_home}/${params.read_glob}"){ file -> file.name.split('_')[0] }
//println("Using read data matching glob: ${params.fastq_home}/*_R{1,2}*.fastq.gz")


//************************
//* Alignment
//************************
process Align {
  tag {"${sample_id}: ${reads[0].getName()} & ${reads[1].getName()}"}
  executor 'local'

  input:
  set sample_id, reads from read_pairs

  output:
  set sample_id, file("${sample_id}-sorted.bam") into aligned_files 

  """
  mkdir -p ${params.alignment_dir}
  mkdir -p ${params.alignment_dir}/logs
  export PERL5LIB
  export PU=`zcat < ${reads[0]} | head -n1 | awk 'BEGIN{FS=":"; OFS="."} /1/{print \$3,\$4,\$10}'`
  
  touch ${sample_id}-sorted.bam
  """
}
// Group files by their sample ids 
joint_bams = aligned_files.groupTuple()


process MergeAligned {
  tag {"${sample_id}"}
  executor 'local'

  input:
  set sample_id, bams from joint_bams

  output:
  set sample_id, file("${sample_id}-merged.bam"), file("${sample_id}-merged.bai") into merged_aligned

  """
  touch ${sample_id}-merged.bam
  touch ${sample_id}-merged.bai
  """
}


process MarkDuplicates {
  publishDir "${params.alignment_dir}", mode: "copy"
  tag {"${sample_id}"}
  executor 'local'

  input:
  set sample_id, file(sample), file(bai) from merged_aligned

  output:
  set val("${sample_id}"), file("${sample_id}-deduped.bam") into deduped

  """
  mkdir -p ${params.alignment_dir}
  mkdir -p ${params.alignment_dir}/logs
  touch ${sample_id}-deduped.bam
  touch ${params.alignment_dir}/logs/${sample_id}.metrix
  """
}
deduped_tumor = deduped.filter{it[0]-params.normal_regex}


// Merge sampes from the same patient into a single bam
process MergeSamples {
  tag {"${sample_id}"}
  executor 'local'
  echo true

  input:
  set sample_id, bam from deduped_tumor

  output:
  set sample_id, file("${sample_id}-merged.bam") into merged_samples

  """
  echo \"Merging: ${params.alignment_dir}/${normal_lookup[sample_id]}-deduped.bam and ${sample_id}-merged.bam\"
  mkdir -p ${params.alignment_dir}
  touch ${sample_id}-merged.bam
  """
}


//************************
//* Co-cleaning 
//************************
// RealignerTargetCreator
// IndelRealigner
// BaseRecalibrator
// PrintReads

// Realign indels
process Indels {
  publishDir "${params.cleaned_dir}", mode: "copy", pattern: "*.intervals"
  tag {"${sample_id}"}
  clusterOptions "vmem=48G,mem=48G"
  executor 'local'

  input:
  set sample_id, file("${sample_id}-merged.bam") from merged_samples

  output:
  set sample_id, file("${sample_id}-indels.bam"), file("${sample_id}-indels.bai") into indeld
  
  """
  mkdir -p ${params.cleaned_dir}
  touch ${sample_id}.intervals
  touch ${sample_id}-indels.bam
  touch ${sample_id}-indels.bai
  """
}


// Recalibrate base quality scores 
process Recal {
  publishDir "${params.cleaned_dir}", mode: "copy", pattern: "*.grp *-clean.bam"
  tag {"${sample_id}"}
  clusterOptions "ncpus=8"
  executor 'local'

  input:
  set sample_id, file("${sample_id}-indels.bam"), file("${sample_id}-indels.bai") from indeld

  output:
  set sample_id, file("${sample_id}-indels.bam"), file("${sample_id}-indels.bai"), file("${sample_id}-bqsr.grp"), file("${normal_lookup[sample_id]}-clean.bam") into recald

  """
  mkdir -p ${params.cleaned_dir}
  touch ${sample_id}-bqsr.grp
  touch ${normal_lookup[sample_id]}-clean.bam
  touch ${normal_lookup[sample_id]}-clean.bai
  """
}


//************************
//* Variant calling 
//************************
process CallVariants {
  publishDir "${params.cleaned_dir}", mode: "copy", pattern: "${sample_id}-clean.bam *.vcf"
  clusterOptions "vmem=64G,mem=64G,ncpus=8"
  tag {"${sample_id}"}
  executor 'local'
  echo true

  input:
  set sample_id, file("${sample_id}-indels.bam"), file("${sample_id}-indels.bai"), file("${sample_id}-bqsr.grp"), file("${normal_lookup[sample_id]}-clean.bam") from recald

  output:
  set sample_id, file("${sample_id}-variants.vcf")

  """
  mkdir -p ${params.variant_dir}
  echo \"Variant detecting with Normal: ${normal_lookup[sample_id]}-clean.bam and Tumor: ${sample_id}-clean.bam\"
  touch ${sample_id}-clean.bam
  touch ${sample_id}-variants.vcf
  """
}


// Given a list of samples and regexs to find patient ids from normal sample names, make a {patient -> [normal:foci]} data structure
def match(samples) {
  lookup = [:]
  // Works like this:
  //  1. For each normal in the flattened sample list a.grep(params.normal_regex).collect{}
  //  2.   Get the patient: a.grep(~/${n-sam}.*/)
  //  3.      And remove the normal sample: -n
  //  4.         And for each of these, map the sample to the [normal, foci] pair .collect{ [(it-params.patient_regex): [n, it]] }
  //  5.      Flattening the nested list sample -> pair list of lists
  //  6.   And merging sample -> pair list maps into a single map
  samples.grep(params.normal_regex).each{ 
    normal -> (samples.grep(~/${normal-params.patient_regex}.*/)-normal).each{ 
      tumor -> lookup[tumor] = normal
    }
  }
  return lookup
}
