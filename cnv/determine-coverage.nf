#!/usr/bin/env nextflow
// Prerequisite: Install cnvkit
// conda config --add channels defaults
// conda config --add channels conda-forge
// conda config --add channels bioconda
// conda install cnvkit

bams1 = Channel.create()
bams2 = Channel.create()
Channel.fromPath("${params.home}/*.bam").separate( bams1, bams2 ){ a -> [a, a] }

/************************
 * Create regions
 ************************/
process MakeAccessibleRegions {
  tag {new File(params.hg19).name}

  output:
  file("accessible.bed") into accessible

  """
  cnvkit.py access ${params.hg19} -x ${params.excludes} -o accessible.bed
  #touch accessible.bed
  """
}

//
process CreateTargets {
  tag {new File(params.panel).name}

  output:
  file("targets.bed") into (targets1, targets2)

  """
  cnvkit.py target ${params.panel} --annotate ${params.refFlat} --split -o targets.bed
  #touch targets.bed
  """
}

//
process CreateAntiTargets {
  tag {new File(params.panel).name}

  input:
  file("targets.bed") from targets1
  file("accessible.bed") from accessible

  output:
  file("antitargets.bed") into antitargets

  """
  cnvkit.py antitarget targets.bed --access accessible.bed -o antitargets.bed
  #touch antitargets.bed
  """
}

//
process ComputeTargetCoverage {
  tag {sample.name.split("-")[0]}
  publishDir "${params.output_dir}", mode: "copy"

  input:
  file("targets.bed") from targets2
  each file(sample) from bams1

  output:
  file("${sample.name.split('-')[0]}.targetcoverage.cnn")

  """
  cnvkit.py coverage ${sample} targets.bed -o ${sample.name.split("-")[0]}.targetcoverage.cnn
  #touch ${sample.name.split("-")[0]}.targetcoverage.cnn
  """
}

//
process ComputeAntiTargetCoverage {
  tag {sample.name.split("-")[0]}
  publishDir "${params.output_dir}", mode: "copy"

  input:
  file("antitargets.bed") from antitargets
  each file(sample) from bams2

  output:
  file("${sample.name.split('-')[0]}.antitargetcoverage.cnn")

  """
  cnvkit.py coverage ${sample} antitargets.bed -o ${sample.name.split("-")[0]}.antitargetcoverage.cnn
  #touch ${sample.name.split("-")[0]}.antitargetcoverage.cnn
  """
}

