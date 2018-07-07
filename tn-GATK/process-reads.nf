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
  clusterOptions "-l vmem=32gb,mem=32gb,nodes=1:ppn=4"
  
  input:
  set sample_id, reads from read_pairs

  output:
  set sample_id, file("${sample_id}-sorted.bam") into aligned_files 

  """
  mkdir -p ${params.alignment_dir}
  mkdir -p ${params.alignment_dir}/logs
  export PERL5LIB
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
  tag {"${sample_id}"}
  clusterOptions "-l vmem=16gb,mem=16gb,nodes=1:ppn=4"

  input:
  set sample_id, bams from joint_bams

  output:
  set sample_id, file("${sample_id}-merged.bam"), file("${sample_id}-merged.bai") into merged_aligned

  """
  java -Xms4g -XX:ParallelGCThreads=4 -jar ${params.picard_home}/picard.jar MergeSamFiles \
    ASSUME_SORTED=true \
    CREATE_INDEX=true \
    INPUT=${bams.join(" INPUT=")} \
    OUTPUT=${sample_id}-merged.bam \
    USE_THREADING=true \
    VALIDATION_STRINGENCY=STRICT
  """
}


process MarkDuplicates {
  publishDir "${params.alignment_dir}", mode: "copy"
  tag {"${sample_id}"}
  clusterOptions "-l vmem=16gb,mem=16gb,nodes=1:ppn=4"

  input:
  set sample_id, file(sample), file(bai) from merged_aligned

  output:
  set sample_id, file("${sample_id}-deduped.bam"), file("${sample_id}-deduped.bai") into deduped

  """
  mkdir -p ${params.alignment_dir}
  mkdir -p ${params.alignment_dir}/logs
  java -Xms4g -XX:ParallelGCThreads=4 -jar ${params.picard_home}/picard.jar MarkDuplicates \
    CREATE_INDEX=true \
    INPUT=${sample} \
    OUTPUT=${sample_id}-deduped.bam \
    M=${params.alignment_dir}/logs/${sample_id}.metrix \
    VALIDATION_STRINGENCY=STRICT \
    REMOVE_DUPLICATES=true
  """
}
deduped_tumor = deduped.filter{it[0]-params.normal_regex}


// Merge tumor/normal sample pairs from the same patient into a single bam
process MergeSamples {
  tag {"${sample_id}"}
  clusterOptions "-l vmem=16gb,mem=16gb,nodes=1:ppn=4"

  input:
  set sample_id, file(sample), file(bai) from deduped_tumor

  output:
  set sample_id, file("${sample_id}-merged.bam"), file("${sample_id}-merged.bai") into merged_samples

  """
  java -Xms4g -XX:ParallelGCThreads=4 -jar ${params.picard_home}/picard.jar MergeSamFiles \
    ASSUME_SORTED=true \
    CREATE_INDEX=true \
    INPUT=${sample} \
    INPUT=${params.alignment_dir}/${normal_lookup[sample_id]}-deduped.bam \
    OUTPUT=${sample_id}-merged.bam \
    USE_THREADING=true \
    VALIDATION_STRINGENCY=STRICT
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
  clusterOptions "-l vmem=48gb,mem=48gb,nodes=1:ppn=4"

  input:
  set sample_id, file(sample), file(bai) from merged_samples

  output:
  set sample_id, file("${sample_id}-indels.bam"), file("${sample_id}-indels.bai") into indeld
  
  """
  mkdir -p ${params.cleaned_dir}
  java -Xms4g -XX:ParallelGCThreads=4 -jar ${params.gatk_home}/GenomeAnalysisTK.jar \
    -T RealignerTargetCreator \
    -R ${params.hg19_reference} \
    -I ${sample} \
    -known ${params.known_indels} \
    --filter_reads_with_N_cigar \
    -nt 24 \
    -o ${sample_id}.intervals

  java -Xms4g -XX:ParallelGCThreads=4 -jar ${params.gatk_home}/GenomeAnalysisTK.jar \
    -T IndelRealigner \
    -R ${params.hg19_reference} \
    -I ${sample} \
    -known ${params.known_indels} \
    -targetIntervals ${sample_id}.intervals \
    --filter_reads_with_N_cigar \
    -o ${sample_id}-indels.bam
  """
}


// Recalibrate base quality scores 
process Recal {
  publishDir "${params.cleaned_dir}", mode: "copy", pattern: "*.grp"
  tag {"${sample_id}"}
  clusterOptions "-l vmem=16gb,mem=16gb,nodes=1:ppn=8"

  input:
  set sample_id, file(sample), file(bai) from indeld

  output:
  set sample_id, file(sample), file(bai), file("${sample_id}-bqsr.grp") into recald

  """
  mkdir -p ${params.cleaned_dir}
  java -Xms4g -XX:ParallelGCThreads=4 -jar ${params.gatk_home}/GenomeAnalysisTK.jar \
    -T BaseRecalibrator \
    -R ${params.hg19_reference} \
    -I ${sample} \
    -knownSites ${params.dbsnps} \
    -nct 8 \
    --filter_reads_with_N_cigar \
    -rf BadCigar \
    -rf MappingQualityZero \
    -o ${sample_id}-bqsr.grp
  """
}


//************************
//* Variant calling 
//************************
process CallVariants {
  publishDir "${params.cleaned_dir}", mode: "copy", pattern: "${sample_id}*"
  clusterOptions "-l vmem=64gb,mem=64gb,nodes=1:ppn=8"
  tag {"${sample_id}"}

  input:
  set sample_id, file(sample), file(bai), file(bqsr) from recald

  output:
  set sample_id, file("${sample_id}.vcf")

  """
  mkdir -p ${params.variant_dir}

  # Make a normal bam with recalibrated bases
  java -Xms4g -XX:ParallelGCThreads=4 -jar ${params.gatk_home}/GenomeAnalysisTK.jar \
    -T PrintReads \
    -R ${params.hg19_reference} \
    -I ${sample} \
    -BQSR ${sample_id}-bqsr.grp \
    -sn ${normal_lookup[sample_id]} \
    -nct 8 \
    -o ${normal_lookup[sample_id]}-clean.bam

  # Make a tumor bam with recalibrated bases
  java -Xms4g -XX:ParallelGCThreads=4 -jar ${params.gatk_home}/GenomeAnalysisTK.jar \
    -T PrintReads \
    -R ${params.hg19_reference} \
    -I ${sample} \
    -BQSR ${sample_id}-bqsr.grp \
    -sn ${sample_id} \
    -nct 8 \
    -o ${sample_id}-clean.bam

  # Call variants 
  /opt/java/jdk1.7.0_latest/bin/java -Xms4g -XX:ParallelGCThreads=4 -jar params.mutect_home/mutect-1.1.7.jar \
     -T MuTect \
     -R ${params.hg19_reference} \
     -I:tumor ${sample_id}-clean.bam \
     -I:normal ${normal_lookup[sample_id]}-clean.bam \
     -o ${sample_id}.tsv
     -vcf ${sample_id}.vcf
  """
}


// Get [tumor: normal] lookup dictionary
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
