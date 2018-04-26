#!/usr/bin/env nextflow

deduped = Channel.fromPath("${params.alignment_dir}/*deduped.bam").map{a -> [ a.name.split("-")[0], a ] }

process CallVariants {
  publishDir "${params.variant_dir}", mode: "copy"
  clusterOptions 'mem=64G,vmem=64G'
  tag {"${sample_id}"}

  input:
  set sample_id, file(sample) from deduped

  output:
  set sample_id, file("${sample_id}-called.vcf")

  """
  mkdir -p ${params.variant_dir}
  ${params.freebayes_home}/freebayes \
    -f ${params.hg19_reference} \
    -F ${params.min_alternate_fraction} \
    -C ${params.min_alternate_count} \
    ${sample} > ${sample_id}-freebayes-called.vcf
  #samtools mpileup -uf ${params.common_variants} ${sample} | bcftools view -Ov > ${sample_id}-called.vcf  #### Check the variant file here
  """
}