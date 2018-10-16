#!/usr/bin/env nextflow
// Prerequisite: Build the genome index like so:
// /wittelab/software/hisat2-2.0.5/hisat2-build /wittelab/hg19/ucsc.hg19.fasta /wittelab/hg19/ucsc.hg19

import groovy.io.FileType

// Get file names (not using a channel)
def fastqs = []
def dir = new File("${params.fastq_home}")
dir.eachFileRecurse (FileType.FILES) { file -> fastqs << file.name.split('_')[0] }
// mapping is a dictionary mapping patient ids to normal, foci pairs
def mapping = match(fastqs.unique().collect())
//println fastqs


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


// Merge cfDNA sample with its normal
// REQUIRES NORMALS PROCESSED TO THIS POINT
process MergeSamples {
  publishDir "${params.alignment_dir}", mode: "copy"
  tag {"${patient_id}"}

  input:
  set patient_id, file("deduped.bam") from deduped

  output:
  set patient_id, file("${patient_id}-merged.bam") into merged_samples

  """
  mkdir -p ${params.alignment_dir}
  java -Xms4g -XX:ParallelGCThreads=4 -jar ${params.picard_home}/picard.jar MergeSamFiles \
    ASSUME_SORTED=true \
    CREATE_INDEX=true \
    INPUT=deduped.bam \
    INPUT= ${find_normal(patient_id)} \
    MERGE_SEQUENCE_DICTIONARIES=false \
    OUTPUT=${patient_id}-merged.bam \
    SORT_ORDER=coordinate \
    USE_THREADING=true \
    VALIDATION_STRINGENCY=STRICT
  """
}


/************************
 * Co-cleaning workflow
 ************************/
// RealignerTargetCreator
// IndelRealigner
// BaseRecalibrator
// PrintReads

// Realign indels
process Indels {
  publishDir "${params.cleaned_dir}", mode: "copy"
  tag {"${patient_id}"}

  input:
  set patient_id, file("${patient_id}-deduped.bam") from merged_samples

  output:
  set patient_id, file("${patient_id}-indels.bam") into indeld
  
  """
  mkdir -p ${params.cleaned_dir}
  mkdir -p ${params.cleaned_dir}/reports
  java -Xms4g -XX:ParallelGCThreads=4 -jar ${params.gatk_home}/GenomeAnalysisTK.jar \
    -T RealignerTargetCreator \
    -R ${params.hg19_reference} \
    -I ${patient_id}-deduped.bam \
    -known ${params.known_indels} \
    --filter_reads_with_N_cigar \
    -o ${patient_id}.intervals

  java -Xms4g -XX:ParallelGCThreads=4 -jar ${params.gatk_home}/GenomeAnalysisTK.jar \
    -T IndelRealigner \
    -R ${params.hg19_reference} \
    -I ${patient_id}-deduped.bam \
    -known ${params.known_indels} \
    -targetIntervals ${patient_id}.intervals \
    --filter_reads_with_N_cigar \
    -o ${patient_id}-indels.bam
    

  cp ${patient_id}.intervals ${params.cleaned_dir}/reports
  """
}


// Recalibrate base quality scores 
process Recal {
  publishDir "${params.cleaned_dir}", mode: "copy"
  tag {"${patient_id}"}

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
    -o ${patient_id}-bqsr.grp

  java -Xms4g -XX:ParallelGCThreads=4 -jar ${params.gatk_home}/GenomeAnalysisTK.jar \
    -T PrintReads \
    -R ${params.hg19_reference} \
    -I ${patient_id}-indels.bam \
    -BQSR ${patient_id}-bqsr.grp \
    -o ${patient_id}-cleaned.bam

  cp ${patient_id}-bqsr.grp ${params.cleaned_dir}/reports
  """
}


/************************
 * Variant calling 
 ************************/
process CallVariants {
  publishDir "${params.variant_dir}", mode: "copy"
  memory '128 GB'
  clusterOptions 'vmem=128 GB'
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
    -o tumor.bam

  # Separate out the normal reads
  java -Xms4g -XX:ParallelGCThreads=4 -jar ${params.gatk_home}/GenomeAnalysisTK.jar \
    -T PrintReads \
    -R ${params.hg19_reference} \
    -sn ${mapping[patient_id][0]} \
    -I ${patient_id}-cleaned.bam \
    -o normal.bam

  # Call variants 
  java -Xms4g -XX:ParallelGCThreads=4 -jar ${params.gatk_home}/GenomeAnalysisTK.jar \
     -T MuTect2 \
     -R ${params.hg19_reference} \
     -I:tumor tumor.bam \
     -I:normal normal.bam \
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


// Given a list of samples and regexs to find patient and normal samples, make a {patient -> [normal:foci]} data structure
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


def find_normal(String id) {
  new File("${params.normal_dir}").eachDirRecurse { dir ->
        dir.eachFileMatch(params.normal_regex) { myfile ->
                return myfile
        }
  }
}