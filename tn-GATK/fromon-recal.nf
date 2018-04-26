#!/usr/bin/env nextflow
// Prerequisite: Build the genome index like so:
// /wittelab/software/hisat2-2.0.5/hisat2-build /wittelab/hg19/ucsc.hg19.fasta /wittelab/hg19/ucsc.hg19

import groovy.io.FileType

// Get file names and put into the fastqs list (not using a nextflow channel)
def fastqs = []
def dir = new File("${params.fastq_home}")
dir.eachFileRecurse (FileType.FILES) { file -> fastqs << file.name.split('_')[0] }
// mapping is a dictionary mapping patient ids to normal, foci pairs
def mapping = match(fastqs.unique().collect())
//println fastqs


// Define the input channel of all read paired samples
indeld = Channel.fromPath("${params.cleaned_dir}"+"/*indels.bam").map({[it.name.split("-")[0], it]})


// Recalibrate base quality scores 
process Recal {
  publishDir "${params.cleaned_dir}", mode: "copy"
  tag {"${patient_id}"}
  clusterOptions "ncpus=8"

  input:
  set patient_id, file("${patient_id}-indels.bam") from indeld

  output:
  set patient_id, file("${patient_id}-cleaned.bam") into recald
   

  """
  mkdir -p ${params.cleaned_dir}
  mkdir -p ${params.cleaned_dir}/reports
  java -Xms4g -XX:ParallelGCThreads=4 -jar ${params.gatk_home}/GenomeAnalysisTK.jar \
    -T BaseRecalibrator \
    -R ${params.hg19_reference} \
    -I ${patient_id}-indels.bam \
    -knownSites ${params.dbsnps} \
    -nct 8 \
    -o ${patient_id}-bqsr.grp

  java -Xms4g -XX:ParallelGCThreads=4 -jar ${params.gatk_home}/GenomeAnalysisTK.jar \
    -T PrintReads \
    -R ${params.hg19_reference} \
    -I ${patient_id}-indels.bam \
    -BQSR ${patient_id}-bqsr.grp \
    -nct 8 \
    -o ${patient_id}-cleaned.bam

  cp ${patient_id}-bqsr.grp ${params.cleaned_dir}/reports
  """
}


/************************
 * Variant calling 
 ***********************/
process CallVariants {
  publishDir "${params.variant_dir}", mode: "copy"
  clusterOptions "vmem=64G,mem=64G,ncpus=8"
  tag {"${patient_id}"}

  input:
  set patient_id, file("${patient_id}-cleaned.bam") from recald

  output:
  set patient_id, file("${patient_id}-called.vcf")

  """
  mkdir -p ${params.variant_dir}

  ## NOTE: samtools likely faster, but would require RG_ID instead of RG_SM
  ## samtools view -bhr RG_ID cleaned.bam > RG_SM.sam

  # Separate out the tumor reads
  java -Xms4g -XX:ParallelGCThreads=4 -jar ${params.gatk_home}/GenomeAnalysisTK.jar \
    -T PrintReads \
    -R ${params.hg19_reference} \
    -sn ${mapping[patient_id][1]} \
    -I ${patient_id}-cleaned.bam \
    -nct 8 \
    -o tumor.bam

  # Separate out the normal reads
  java -Xms4g -XX:ParallelGCThreads=4 -jar ${params.gatk_home}/GenomeAnalysisTK.jar \
    -T PrintReads \
    -R ${params.hg19_reference} \
    -sn ${mapping[patient_id][0]} \
    -I ${patient_id}-cleaned.bam \
    -nct 8 \
    -o normal.bam

  # Call variants 
  java -Xms4g -XX:ParallelGCThreads=4 -jar ${params.gatk_home}/GenomeAnalysisTK.jar \
     -T MuTect2 \
     -R ${params.hg19_reference} \
     -I:tumor tumor.bam \
     -I:normal normal.bam \
     -nct 8 \
     -o ${patient_id}-called.vcf
  """
}


// Get a base filename and optionally add a new extention
def base(file, ext = "") {
  name = file.name
  name = name.take(name.lastIndexOf('.'))
  //name = name.take(name.indexOf('-'))
  name = name+ext
  return name
}


// Given a list of samples and regexs to find patient ids from normal sample names, make a {patient -> [normal:foci]} data structure
def match(a) {
  // Works like this:
  //  1. For each normal in the flattened sample list a.grep(params.normal_regex).collect{}
  //  2.   Get the patient: a.grep(~/${n-sam}.*/)
  //  3.      And remove the normal sample: -n
  //  4.         And for each of these, map the sample to the [normal, foci] pair .collect{ [(it-params.patient_regex): [n, it]] }
  //  5.      Flattening the nested list sample -> pair list of lists
  //  6.   And merging sample -> pair list maps into a single map
  return a.grep(params.normal_regex).collect{ 
    n -> (a.grep(~/${n-params.patient_regex}.*/)-n).collect{ 
      [(it-params.patient_regex): [n, it]] 
    }
  }.flatten().collectEntries()
}
