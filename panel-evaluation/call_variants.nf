#!/usr/bin/env nextflow
deduped = Channel.fromPath("${params.dedup_dir}/P*.bam").map{file -> [file.name.split('-')[0], file, file+".bai"] }

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

/* Other methods

process CallVariants {
  echo true
  beforeScript 'source /home/carioc/.zshrc_local'
  publishDir "${params.variant_dir}", mode: "copy", pattern: "*.vcf"
  tag {"${sample_id}"}
  memory '256 GB'

  input:
  set sample_id, file(sample), file(sample_index) from deduped

  output:
  set sample_id, file("${sample_id}-called.vcf")

  """
  mkdir -p ${params.variant_dir}
  samtools index ${sample}
  freebayes-1.0 -f ${params.hg19_reference} ${sample} > ${sample_id}-called.vcf
  #java -jar /wittelab/software/GATK/picard.jar AddOrReplaceReadGroups \
  # I=${sample} \
  # O=rg-${sample} \
  # LB=cfDNA \
  # PL=illumina \
  # PU=cfDNA \
  # SM=cfDNA
  #samtools index rg-${sample}
  #java -jar /wittelab/software/GATK/GenomeAnalysisTK.jar \
  # -T UnifiedGenotyper \
  # -R ${params.hg19_reference} \
  # -glm BOTH \
  # -I rg-${sample} \
  # -o ${sample_id}-called.vcf \
  # -U ALLOW_N_CIGAR_READS
  # #[-L targets.interval_list]
  """
}

*/
