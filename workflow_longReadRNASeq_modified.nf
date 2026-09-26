#!/usr/bin/env nextflow
// Based on: https://github.com/GoekeLab/sg-nex-data/blob/master/docs/colab/workflow_longReadRNASeq.nf
// Modified for DSA4262 Assignment 1 Task 5:
//   - protocol-specific minimap2 flags (direct RNA vs cDNA), read from a samplesheet
//   - optional genome annotations for Bambu (params.refGtf may be omitted -> discovery without annotations)
//   - added a QC process (samtools flagstat-based mapping rate check)
//   - Bambu is run jointly on all samples (collected bam files) instead of one bam at a time

params.samplesheet = '/path/to/samples.csv'   // columns: sample_id,fastq_path,protocol
params.refFa       = '/path/to/ref.fa'
params.refGtf      = null                     // omit to run Bambu WITHOUT annotations
params.outdir      = 'results'
params.qc_min_mapping_rate = 50               // percent

process MINIMAP2_ALIGN {
  tag "$sample_id"
  cpus 6

  input:
    tuple val(sample_id), path(reads), val(protocol)
    path refFa

  output:
    tuple val(sample_id), path("${sample_id}.sam")

  script:
  def protocol_flags = (protocol == 'directRNA') ? '-uf -k14' : ''
  """
  minimap2 -ax splice ${protocol_flags} -t ${task.cpus} ${refFa} ${reads} > ${sample_id}.sam
  """
}

process SAM_TO_BAM {
  tag "$sample_id"
  cpus 4
  publishDir "${params.outdir}/bam", mode: 'copy'

  input:
    tuple val(sample_id), path(sam)

  output:
    tuple val(sample_id), path("${sample_id}.bam"), path("${sample_id}.bam.bai")

  script:
  """
  samtools sort -@ ${task.cpus} -o ${sample_id}.bam ${sam}
  samtools index ${sample_id}.bam
  """
}

process QC {
  tag "$sample_id"
  publishDir "${params.outdir}/qc", mode: 'copy'

  input:
    tuple val(sample_id), path(bam), path(bai)

  output:
    path "${sample_id}.qc_report.txt"

  script:
  """
  samtools flagstat ${bam} > ${sample_id}.flagstat.txt
  total=\$(samtools view -c ${bam})
  mapped=\$(samtools view -c -F 4 ${bam})
  rate=\$(awk -v m="\$mapped" -v t="\$total" 'BEGIN{ if (t==0) print 0; else printf "%.2f", (m/t)*100 }')
  status=\$(awk -v r="\$rate" -v thr="${params.qc_min_mapping_rate}" 'BEGIN{ print (r>=thr) ? "PASS" : "FAIL" }')
  {
    echo "sample: ${sample_id}"
    echo "total_reads: \$total"
    echo "mapped_reads: \$mapped"
    echo "mapping_rate_pct: \$rate"
    echo "qc_threshold_pct: ${params.qc_min_mapping_rate}"
    echo "QC_STATUS: \$status"
  } > ${sample_id}.qc_report.txt
  cat ${sample_id}.flagstat.txt >> ${sample_id}.qc_report.txt
  """
}

process BAMBU {
  cpus 6
  publishDir params.outdir, mode: 'copy'

  input:
    path refFa
    path bams
    path bais
    val refGtfPath

  output:
    path "counts_transcript.txt"
    path "counts_gene.txt"
    path "extended_annotations.gtf"

  script:
  def annotation_code = (refGtfPath != 'NONE') ?
      "annotations <- prepareAnnotations('${refGtfPath}')" :
      "annotations <- NULL"
  """
  #!/home/ubuntu/miniforge3/envs/bambu/bin/Rscript --vanilla
  library(bambu)
  ${annotation_code}
  bam_files <- strsplit("${bams}", " ")[[1]]
  se <- bambu(reads = bam_files, annotations = annotations, genome = "${refFa}", ncore = ${task.cpus})
  writeBambuOutput(se, path = "./")
  """
}

workflow {
  samples_ch = Channel
    .fromPath(params.samplesheet)
    .splitCsv(header: true)
    .map { row -> tuple(row.sample_id, file(row.fastq_path), row.protocol) }

  aligned_ch = MINIMAP2_ALIGN(samples_ch, params.refFa)
  bam_ch     = SAM_TO_BAM(aligned_ch)

  QC(bam_ch)

  bam_paths = bam_ch.map { it[1] }.collect()
  bai_paths = bam_ch.map { it[2] }.collect()

  BAMBU(params.refFa, bam_paths, bai_paths, params.refGtf ?: 'NONE')
}
