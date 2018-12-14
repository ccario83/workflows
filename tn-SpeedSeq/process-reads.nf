#!/usr/bin/env nextflow
read_pairs = Channel.fromFilePairs("${params.fastq_home}/${params.read_glob}"){ file -> file.name.split('_')[0] }

//************************
//* Alignment
//************************
process SpeedAlign {
  tag {"${sample_id}: ${reads[0].getName()} & ${reads[1].getName()}"}
  clusterOptions "-l vmem=64gb,mem=64gb,nodes=1:ppn=12"

  input:
  set sample_id, reads from read_pairs

  output:
  set sample_id, file("${sample_id}.bam") into aligned_files 

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
// Group files by their sample ids 
joint_bams = aligned_files.groupTuple()

process SpeedMerge {
  publishDir "${params.alignment_dir}", mode: "copy", pattern: "*.bam *.bai"
  tag {"${sample_id}"}
  clusterOptions "-l vmem=32gb,mem=32gb,nodes=1:ppn=12"

  input:
  set sample_id, bams from joint_bams

  output:
  set sample_id, file("${sample_id}_merged.bam"), file("${sample_id}_merged.bam.bai") into merged_aligned

  script:
  if(bams.size()==1)
    """
    mkdir -p ${params.alignment_dir}

    ln -s ${bams[0]} ${sample_id}_merged.bam 
    $SAMBAMBA index -t 8 ${sample_id}_merged.bam
    """
  else
    """
    mkdir -p ${params.alignment_dir}

    $SAMBAMBA merge -t 10 ${sample_id}_merged.bam ${bams.each{it[0]}.join(" ")}
    $SAMBAMBA index -t 10 ${sample_id}_merged.bam
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
process SpeedCall {
  publishDir "${params.variant_dir}", mode: "copy", pattern: "*.vcf.gz"
  tag {"${tumor[0]}"}
  clusterOptions "-l vmem=100gb,mem=100gb,nodes=1:ppn=32"


  input:
  set normal, tumor from paired

  output:
  file("${tumor[0]}.vcf.gz")

  """
  mkdir -p ${params.variant_dir}
  export LD_LIBRARY_PATH="/opt/gcc/gcc-4.8.1/lib64/"
  
  $SPEEDSEQ somatic \
    -o ${tumor[0]} \
    -F 0.05 \
    -q 1 \
    -t 30 \
    ${params.hg19_reference} \
    ${normal[1]} \
    ${tumor[1]}
  """
}

