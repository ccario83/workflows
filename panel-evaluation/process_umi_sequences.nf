#!/usr/bin/env nextflow
// Prerequisite: Build the genome index like so:
// /wittelab/software/hisat2-2.0.5/hisat2-build /wittelab/hg19/ucsc.hg19.fasta /wittelab/hg19/ucsc.hg19

// Define the input channel
//Channel.fromFilePairs("${params.fastq_home}/*_R{1,2}*.fastq"){ file -> file.name.split('_')[0] }.println { id, files -> "Files with the sample ID $id are $files" }
read_pairs = Channel.fromFilePairs("${params.fastq_home}/*_R{1,2}*.fastq"){ file -> file.name.split('_')[0] }
println("${params.fastq_home}/*_R{1,2}*.fastq")

process Align {
    tag {"${sample_id}: ${reads[0]} & ${reads[1]}"}
    cache true

    input:
    set sample_id, reads from read_pairs
 
    output:
    set sample_id, file("${sample_id}-sorted.bam") into aligned_files 
 
    """
    mkdir -p ${params.alignment_dir}
    mkdir -p ${params.alignment_dir}/logs
    ${params.hisat2_home}/hisat2 -f -x ${params.hisat2_index} -q -1 ${reads[0]} -2 ${reads[1]} -S ${sample_id}.sam &> ${params.alignment_dir}/logs/${base(reads[0], ".log")}
    samtools view -bS ${sample_id}.sam | samtools sort -@ 1 -o ${sample_id}-sorted.bam
    """
}
joint_bams = aligned_files.groupTuple()



process MergeAligned {
  publishDir "${params.alignment_dir}", mode: "copy"
  tag {"${sample_id}"}
  cache true

  input:
  set sample_id, bams from joint_bams

  output:
  set sample_id, file("${sample_id}-merged.bam"), file("${sample_id}-merged.bam.bai") into merged
   

  """
  mkdir -p ${params.alignment_dir}
  samtools merge -f ${sample_id}-merged.bam ${bams.join(" ")}
  samtools index ${sample_id}-merged.bam
  """
}

/* Base recalibration
process Recalibrate {
  publishDir "${params.recal_dir}", mode: "copy"

  tag {"${sample_id}"}

  input:
  set sample_id, file(sample), file(sample_index) from merged

  output:
  set sample_id, file("${sample_id}-recalibrated.bam"), file("${sample_id}-recalibrated.bam.bai") into recald
   

  """
  mkdir -p ${params.recal_dir}
  java -jar GenomeAnalysisTK.jar \
     -T PrintReads \
     -R params.hg19_reference \
     -I input.bam \
     -BQSR recalibration_report.grp \
     -o ${sample_id}-recalibrated.bam
  samtools index ${sample_id}-recalibrated.bam
  """

}
*/

process ProcessUMIs {
  errorStrategy 'ignore'
  publishDir "${params.dedup_dir}", mode: "copy"
  tag {"${sample_id}"}

  input:
  set sample_id, file(sample), file(sample_index) from merged

  output:
  set sample_id, file("${sample_id}-deduped.bam"), file("${sample_id}-deduped.bam.bai")  into deduped

  """
  mkdir -p ${params.dedup_dir}
  connor -f ${params.consensus_freq_threshold}\
  -s ${params.min_family_size_threshold}\
  -d {params.umt_distance_threshold}\
  ${sample} ${sample_id}-deduped.bam --force
  """
}

process CallVariants {
  beforeScript 'source /home/carioc/.zshrc_local'
  publishDir "${params.variant_dir}", mode: "copy"
  memory '256 GB'
  tag {"${sample_id}"}

  input:
  set sample_id, file(sample), file(sample_index) from deduped

  output:
  set sample_id, file("${sample_id}-called.vcf")

  """
  mkdir -p ${params.variant_dir}
  #samtools index ${sample}
  freebayes-1.0\
  -f ${params.hg19_reference}\
  -F ${params.min_alternate_fraction}\
  -C ${params.min_alternate_count}\
  ${sample} > ${sample_id}-called.vcf
  #samtools mpileup -uf ${params.common_variants} ${sample} | bcftools view -Ov > ${sample_id}-called.vcf  #### Check the variant file here
  """
}

def base(file, ext = "") {
  name = file.name
  name = name.take(name.lastIndexOf('.'))
  //name = name.take(name.indexOf('-'))
  name = name+ext
  return name
}
