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
<<<<<<< HEAD
    -d ${params.umt_distance_threshold}\
=======
    -d {params.umt_distance_threshold}\
>>>>>>> 05fbfb22bc07e6bbd9ecfefd151e4100c73c952d
    --force \
    ${sample} \
    ${sample_id}-collapsed.bam 
  """
}


//************************
<<<<<<< HEAD
//* Variant calling with Freebayes or VarDict
//************************
process Call {
  tag {"${sample_id} w/ ${params.variant_caller}"}
  clusterOptions "-l vmem=64gb,mem=64gb"
  publishDir "${params.variant_dir}", mode: "copy", pattern: "*.vcf*"
=======
//* Variant calling 
//************************
process Call {
  tag {"${sample_id}"}
  clusterOptions "-l vmem=256gb,mem=256gb"
  publishDir "${params.variant_dir}", mode: "copy", pattern: "*.vcf.gz"
>>>>>>> 05fbfb22bc07e6bbd9ecfefd151e4100c73c952d
  

  input:
  set sample_id, file(sample), file(sample_index) from collapsed

  output:
  set sample_id, file("${sample_id}.vcf")

<<<<<<< HEAD
  script:
  if (params.variant_caller == 'freebayes')
    """
    mkdir -p ${params.variant_dir}
    
    ${params.freebayes_home}/freebayes \
      -f ${params.hg19_reference}\
      --min-alternate-fraction ${params.min_alternate_fraction}\
      --min-alternate-count ${params.min_alternate_count}\
      ${sample} > ${sample_id}.vcf
    #cat ${sample_id}.vcf | ${params.vcffilter_home}/vcffilter -f "${params.variant_filter}" > ${sample_id}-filtered.vcf
    """
  else if (params.variant_caller == 'vardict')
    """
    mkdir -p ${params.variant_dir}
    
    ${params.vardict_home}/VarDict \
      -G ${params.hg19_reference}\
      -f ${params.min_alternate_fraction}\
      -N ${sample_id}\
      -b ${sample}\
      -z 0\
      -c 1\
      -S 2\
      -E 3\
      -g 4\
      ${params.panel}\
      | ${params.vardict_script_home}/teststrandbias.R \
      | ${params.vardict_script_home}/var2vcf_valid.pl -N ${sample_id} -E -f ${params.min_alternate_fraction} \
      > ${sample_id}.vcf
    """
}



=======
  """
  mkdir -p ${params.variant_dir}
  
  freebayes-1.0\
    -f ${params.hg19_reference}\
    -F ${params.min_alternate_fraction}\
    -C ${params.min_alternate_count}\
    ${sample} > ${sample_id}.vcf
  """
}


>>>>>>> 05fbfb22bc07e6bbd9ecfefd151e4100c73c952d
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
