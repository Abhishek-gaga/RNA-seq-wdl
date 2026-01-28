version 1.0

workflow borrelia_rnaseq_umi {
  input {
    # --- Required per-sample inputs ---
    File   demux_bam            # unaligned demux BAM (FASTQ-in-BAM)
    String sample_id

    # --- Reference inputs ---
    String bowtie2_index_prefix # e.g. gs://bucket/index/borrelia (prefix, not .bt2 file)
    File   annotation_gtf       # gene annotation (GTF recommended)
    Int    strandedness = 0     # featureCounts: 0 unstranded, 1 stranded, 2 reversely stranded

    # --- Trimming inputs ---
    File?  trimmomatic_adapters # optional; if null we'll use a safe default adapter file
    Int    minlen = 30
    Int    threads = 8

    # --- Step-by-step control ---
    # options: "bam_to_fastq", "fastqc_raw", "trim", "fastqc_trim", "align",
    #          "tag_umi", "dedup", "count", "multiqc", "none"
    String stop_after = "none"

    # --- Runtime knobs ---
    Int  disk_gb = 200
    Int  memory_gb = 32
  }

  call BamToFastqWithUMI {
    input:
      demux_bam = demux_bam,
      sample_id = sample_id,
      threads   = threads,
      disk_gb   = disk_gb,
      memory_gb = memory_gb
  }

  if (stop_after == "bam_to_fastq") {
    output {
      File out_r1_fastq = BamToFastqWithUMI.r1_fastq_gz
      File out_r2_fastq = BamToFastqWithUMI.r2_fastq_gz
      File umi_rewrite_log = BamToFastqWithUMI.umi_rewrite_log
    }
  }

  call FastQC as FastQC_Raw {
    input:
      sample_id = sample_id + ".raw",
      r1_fastq_gz = BamToFastqWithUMI.r1_fastq_gz,
      r2_fastq_gz = BamToFastqWithUMI.r2_fastq_gz,
      threads = threads,
      disk_gb = disk_gb,
      memory_gb = memory_gb
  }

  if (stop_after == "fastqc_raw") {
    output {
      Array[File] raw_fastqc = FastQC_Raw.fastqc_zip
      Array[File] raw_fastqc_html = FastQC_Raw.fastqc_html
    }
  }

  call TrimmomaticPE {
    input:
      sample_id = sample_id,
      r1_fastq_gz = BamToFastqWithUMI.r1_fastq_gz,
      r2_fastq_gz = BamToFastqWithUMI.r2_fastq_gz,
      adapters = trimmomatic_adapters,
      threads  = threads,
      minlen   = minlen,
      disk_gb  = disk_gb,
      memory_gb= memory_gb
  }

  if (stop_after == "trim") {
    output {
      File trim_r1 = TrimmomaticPE.r1_paired_fastq_gz
      File trim_r2 = TrimmomaticPE.r2_paired_fastq_gz
      File trim_log = TrimmomaticPE.trim_log
    }
  }

  call FastQC as FastQC_Trim {
    input:
      sample_id = sample_id + ".trim",
      r1_fastq_gz = TrimmomaticPE.r1_paired_fastq_gz,
      r2_fastq_gz = TrimmomaticPE.r2_paired_fastq_gz,
      threads = threads,
      disk_gb = disk_gb,
      memory_gb = memory_gb
  }

  if (stop_after == "fastqc_trim") {
    output {
      Array[File] trim_fastqc = FastQC_Trim.fastqc_zip
      Array[File] trim_fastqc_html = FastQC_Trim.fastqc_html
    }
  }

  call Bowtie2Align {
    input:
      sample_id = sample_id,
      bowtie2_index_prefix = bowtie2_index_prefix,
      r1_fastq_gz = TrimmomaticPE.r1_paired_fastq_gz,
      r2_fastq_gz = TrimmomaticPE.r2_paired_fastq_gz,
      threads = threads,
      disk_gb = disk_gb,
      memory_gb = memory_gb
  }

  if (stop_after == "align") {
    output {
      File aligned_bam = Bowtie2Align.sorted_bam
      File aligned_bai = Bowtie2Align.sorted_bai
      File align_log   = Bowtie2Align.align_log
      File flagstat    = Bowtie2Align.flagstat_txt
    }
  }

  # Convert UMI embedded in read_id -> RX tag in BAM (so umi_tools dedup can use it)
  call UmiToolsTagFromReadID {
    input:
      sample_id = sample_id,
      in_bam = Bowtie2Align.sorted_bam,
      disk_gb = disk_gb,
      memory_gb = memory_gb
  }

  if (stop_after == "tag_umi") {
    output {
      File tagged_bam = UmiToolsTagFromReadID.tagged_bam
      File tagged_bai = UmiToolsTagFromReadID.tagged_bai
      File tag_log    = UmiToolsTagFromReadID.tag_log
    }
  }

  call UmiToolsDedup {
    input:
      sample_id = sample_id,
      in_bam = UmiToolsTagFromReadID.tagged_bam,
      disk_gb = disk_gb,
      memory_gb = memory_gb
  }

  if (stop_after == "dedup") {
    output {
      File dedup_bam = UmiToolsDedup.dedup_bam
      File dedup_bai = UmiToolsDedup.dedup_bai
      File dedup_log = UmiToolsDedup.dedup_log
    }
  }

  call FeatureCounts {
    input:
      sample_id = sample_id,
      in_bam = UmiToolsDedup.dedup_bam,
      annotation_gtf = annotation_gtf,
      strandedness = strandedness,
      threads = threads,
      disk_gb = disk_gb,
      memory_gb = memory_gb
  }

  if (stop_after == "count") {
    output {
      File counts_tsv = FeatureCounts.counts_tsv
      File counts_log = FeatureCounts.counts_log
    }
  }

  # MultiQC to summarize FastQC + (optionally) alignment logs
  call MultiQC {
    input:
      sample_id = sample_id,
      fastqc_zip = FastQC_Raw.fastqc_zip + FastQC_Trim.fastqc_zip,
      fastqc_html = FastQC_Raw.fastqc_html + FastQC_Trim.fastqc_html,
      extra_files = [Bowtie2Align.flagstat_txt, Bowtie2Align.align_log, TrimmomaticPE.trim_log, UmiToolsDedup.dedup_log],
      disk_gb = disk_gb,
      memory_gb = memory_gb
  }

  if (stop_after == "multiqc") {
    output {
      File multiqc_html = MultiQC.multiqc_html
      File multiqc_data = MultiQC.multiqc_data_zip
    }
  }

  output {
    # Final “full run” outputs
    File r1_fastq_gz = BamToFastqWithUMI.r1_fastq_gz
    File r2_fastq_gz = BamToFastqWithUMI.r2_fastq_gz

    File trim_r1_fastq_gz = TrimmomaticPE.r1_paired_fastq_gz
    File trim_r2_fastq_gz = TrimmomaticPE.r2_paired_fastq_gz

    File aligned_bam = Bowtie2Align.sorted_bam
    File aligned_bai = Bowtie2Align.sorted_bai
    File dedup_bam   = UmiToolsDedup.dedup_bam
    File dedup_bai   = UmiToolsDedup.dedup_bai

    File counts_tsv  = FeatureCounts.counts_tsv

    Array[File] fastqc_raw_zip  = FastQC_Raw.fastqc_zip
    Array[File] fastqc_trim_zip = FastQC_Trim.fastqc_zip
    File multiqc_html = MultiQC.multiqc_html

    # Logs (useful for production)
    File umi_rewrite_log = BamToFastqWithUMI.umi_rewrite_log
    File trim_log        = TrimmomaticPE.trim_log
    File align_log       = Bowtie2Align.align_log
    File flagstat_txt    = Bowtie2Align.flagstat_txt
    File tag_log         = UmiToolsTagFromReadID.tag_log
    File dedup_log       = UmiToolsDedup.dedup_log
    File counts_log      = FeatureCounts.counts_log
  }
}

# -----------------------------
# Tasks
# -----------------------------

task BamToFastqWithUMI {
  input {
    File   demux_bam
    String sample_id
    Int    threads
    Int    disk_gb
    Int    memory_gb
  }

  command <<<
    set -euo pipefail
    samtools fastq -@ ~{threads} \
      -1 ~{sample_id}_R1.fastq.gz \
      -2 ~{sample_id}_R2.fastq.gz \
      -0 /dev/null -s /dev/null -n \
      umi_in_qname.bam

  >>>

  output {
    File r1_fastq_gz = "~{sample_id}_R1.fastq.gz"
    File r2_fastq_gz = "~{sample_id}_R2.fastq.gz"
    File umi_rewrite_log = "umi_rewrite.log"
  }

  runtime {
    docker: "yourdockerhub/bacterial-rnaseq-umi:1.0.0"
    cpu: threads
    memory: "~{memory_gb} GB"
    disks: "local-disk ~{disk_gb} HDD"
  }
}

task FastQC {
  input {
    String sample_id
    File   r1_fastq_gz
    File   r2_fastq_gz
    Int    threads
    Int    disk_gb
    Int    memory_gb
  }

  command <<<
    set -euo pipefail
    mkdir -p fastqc_out
    fastqc -t ~{threads} -o fastqc_out ~{r1_fastq_gz} ~{r2_fastq_gz}
    ls -lh fastqc_out
  >>>

  output {
    Array[File] fastqc_zip  = glob("fastqc_out/*.zip")
    Array[File] fastqc_html = glob("fastqc_out/*.html")
  }

  runtime {
    docker: "yourdockerhub/bacterial-rnaseq-umi:1.0.0"
    cpu: threads
    memory: "~{memory_gb} GB"
    disks: "local-disk ~{disk_gb} HDD"
  }
}

task TrimmomaticPE {
  input {
    String sample_id
    File   r1_fastq_gz
    File   r2_fastq_gz
    File?  adapters
    Int    threads
    Int    minlen
    Int    disk_gb
    Int    memory_gb
  }

  command <<<
    set -euo pipefail

    # If adapters not provided, create a minimal adapter file.
    # (You can replace this with a proper Illumina adapter file in your bucket.)
    if [ "~{true=defined(adapters) false=''}" = "true" ]; then
      cp "~{adapters}" adapters.fa
    else
      cat > adapters.fa <<'EOF'
>TruSeq3-PE-2
AGATCGGAAGAGCACACGTCTGAACTCCAGTCA
EOF
    fi

    trimmomatic PE -threads ~{threads} -phred33 \
      ~{r1_fastq_gz} ~{r2_fastq_gz} \
      ~{sample_id}.R1.paired.fastq.gz ~{sample_id}.R1.unpaired.fastq.gz \
      ~{sample_id}.R2.paired.fastq.gz ~{sample_id}.R2.unpaired.fastq.gz \
      ILLUMINACLIP:adapters.fa:2:30:10 \
      LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:~{minlen} \
      2> ~{sample_id}.trimmomatic.log

  >>>

  output {
    File r1_paired_fastq_gz = "~{sample_id}.R1.paired.fastq.gz"
    File r2_paired_fastq_gz = "~{sample_id}.R2.paired.fastq.gz"
    File trim_log = "~{sample_id}.trimmomatic.log"
  }

  runtime {
    docker: "yourdockerhub/bacterial-rnaseq-umi:1.0.0"
    cpu: threads
    memory: "~{memory_gb} GB"
    disks: "local-disk ~{disk_gb} HDD"
  }
}

task Bowtie2Align {
  input {
    String sample_id
    String bowtie2_index_prefix
    File   r1_fastq_gz
    File   r2_fastq_gz
    Int    threads
    Int    disk_gb
    Int    memory_gb
  }

  command <<<
    set -euo pipefail

    bowtie2 -x "~{bowtie2_index_prefix}" \
      -1 ~{r1_fastq_gz} -2 ~{r2_fastq_gz} \
      -p ~{threads} \
      2> ~{sample_id}.bowtie2.log \
    | samtools sort -@ ~{threads} -o ~{sample_id}.sorted.bam

    samtools index ~{sample_id}.sorted.bam
    samtools flagstat ~{sample_id}.sorted.bam > ~{sample_id}.flagstat.txt
  >>>

  output {
    File sorted_bam   = "~{sample_id}.sorted.bam"
    File sorted_bai   = "~{sample_id}.sorted.bam.bai"
    File align_log    = "~{sample_id}.bowtie2.log"
    File flagstat_txt = "~{sample_id}.flagstat.txt"
  }

  runtime {
    docker: "yourdockerhub/bacterial-rnaseq-umi:1.0.0"
    cpu: threads
    memory: "~{memory_gb} GB"
    disks: "local-disk ~{disk_gb} HDD"
  }
}

task UmiToolsTagFromReadID {
  input {
    String sample_id
    File   in_bam
    Int    disk_gb
    Int    memory_gb
  }

  command <<<
    set -euo pipefail

    # We embedded UMI as "..._UMI:ACGT..." in the read name.
    # Now tag it back into BAM as RX so umi_tools dedup can use it.
    umi_tools tag \
      -I ~{in_bam} \
      -S ~{sample_id}.tagged.bam \
      --extract-umi-method=read_id \
      --umi-separator="_UMI:" \
      --umi-tag=RX \
      2> ~{sample_id}.umi_tag.log

    samtools index ~{sample_id}.tagged.bam
  >>>

  output {
    File tagged_bam = "~{sample_id}.tagged.bam"
    File tagged_bai = "~{sample_id}.tagged.bam.bai"
    File tag_log    = "~{sample_id}.umi_tag.log"
  }

  runtime {
    docker: "yourdockerhub/bacterial-rnaseq-umi:1.0.0"
    cpu: 2
    memory: "~{memory_gb} GB"
    disks: "local-disk ~{disk_gb} HDD"
  }
}

task UmiToolsDedup {
  input {
    String sample_id
    File   in_bam
    Int    disk_gb
    Int    memory_gb
  }

  command <<<
    set -euo pipefail

    umi_tools dedup \
      -I ~{in_bam} \
      -S ~{sample_id}.dedup.bam \
      --paired \
      --umi-tag=RX \
      2> ~{sample_id}.umi_dedup.log

    samtools index ~{sample_id}.dedup.bam
  >>>

  output {
    File dedup_bam = "~{sample_id}.dedup.bam"
    File dedup_bai = "~{sample_id}.dedup.bam.bai"
    File dedup_log = "~{sample_id}.umi_dedup.log"
  }

  runtime {
    docker: "yourdockerhub/bacterial-rnaseq-umi:1.0.0"
    cpu: 2
    memory: "~{memory_gb} GB"
    disks: "local-disk ~{disk_gb} HDD"
  }
}

task FeatureCounts {
  input {
    String sample_id
    File   in_bam
    File   annotation_gtf
    Int    strandedness
    Int    threads
    Int    disk_gb
    Int    memory_gb
  }

  command <<<
    set -euo pipefail

    # -p for paired-end
    # -s strandedness (0/1/2)
    featureCounts \
      -a ~{annotation_gtf} \
      -o ~{sample_id}.geneCounts.tsv \
      -T ~{threads} \
      -p \
      -s ~{strandedness} \
      ~{in_bam} \
      > ~{sample_id}.featureCounts.log
  >>>

  output {
    File counts_tsv = "~{sample_id}.geneCounts.tsv"
    File counts_log = "~{sample_id}.featureCounts.log"
  }

  runtime {
    docker: "yourdockerhub/bacterial-rnaseq-umi:1.0.0"
    cpu: threads
    memory: "~{memory_gb} GB"
    disks: "local-disk ~{disk_gb} HDD"
  }
}

task MultiQC {
  input {
    String sample_id
    Array[File] fastqc_zip
    Array[File] fastqc_html
    Array[File] extra_files
    Int    disk_gb
    Int    memory_gb
  }

  command <<<
    set -euo pipefail
    mkdir -p qc
    # Put files in one directory so multiqc picks them up reliably
    for f in ~{sep=' ' fastqc_zip} ~{sep=' ' fastqc_html} ~{sep=' ' extra_files}; do
      cp "$f" qc/ || true
    done
    multiqc qc -o qc
  >>>

  output {
    File multiqc_html     = "qc/multiqc_report.html"
    File multiqc_data_zip = "qc/multiqc_data.zip"
  }

  runtime {
    docker: "yourdockerhub/bacterial-rnaseq-umi:1.0.0"
    cpu: 1
    memory: "~{memory_gb} GB"
    disks: "local-disk ~{disk_gb} HDD"
  }
}
