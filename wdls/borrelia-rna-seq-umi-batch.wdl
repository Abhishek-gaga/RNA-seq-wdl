version 1.0
import "borrelia_rnaseq_umi.wdl" as single

workflow borrelia_rnaseq_umi_batch {
  input {
    # TSV header:
    # sample_id  demux_bam  bowtie2_index_prefix  annotation_gtf  strandedness
    File sample_tsv

    # Global controls (shared)
    File? trimmomatic_adapters
    Int   minlen = 30
    Int   threads = 8
    String stop_after = "none"

    # Runtime knobs (global defaults)
    Int disk_gb = 200
    Int memory_gb = 32
  }

  Array[Array[String]] rows = read_tsv(sample_tsv)

  scatter (i in range(1, length(rows))) {
    String sample_id            = rows[i][0]
    File   demux_bam            = rows[i][1]
    String bowtie2_index_prefix = rows[i][2]
    File   annotation_gtf       = rows[i][3]
    Int    strandedness         = rows[i][4]  # needs to be numeric in TSV

    call single.borrelia_rnaseq_umi as per_sample {
      input:
        demux_bam            = demux_bam,
        sample_id            = sample_id,
        bowtie2_index_prefix = bowtie2_index_prefix,
        annotation_gtf       = annotation_gtf,
        strandedness         = strandedness,
        trimmomatic_adapters = trimmomatic_adapters,
        minlen               = minlen,
        threads              = threads,
        stop_after           = stop_after,
        disk_gb              = disk_gb,
        memory_gb            = memory_gb
    }
  }

  output {
    Array[String] sample_ids = per_sample.sample_id
    Array[File]   counts_tsv = per_sample.counts_tsv
    Array[File]   multiqc_html = per_sample.multiqc_html
    Array[File]   aligned_bam = per_sample.aligned_bam
    Array[File]   dedup_bam = per_sample.dedup_bam
    Array[File]   trim_log = per_sample.trim_log
    Array[File]   align_log = per_sample.align_log
    Array[File]   dedup_log = per_sample.dedup_log
    Array[File]   flagstat_txt = per_sample.flagstat_txt
  }
}
