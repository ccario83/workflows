merged = Channel.fromPath("${params.alignment_dir}/P*.bam").map{file -> [file.name.split('-')[0], file, file+".bai"] }


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
