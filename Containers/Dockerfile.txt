# Bacterial RNA-seq (+UMI) toolbox container
# Includes: fastqc, multiqc, trimmomatic, bowtie2, samtools, umi_tools, featureCounts (subread)

FROM debian:bookworm-slim

# Optional: use uv for fast Python venv + installs (like your earlier Dockerfile)
COPY --from=ghcr.io/astral-sh/uv:0.7.22 /uv /uvx /bin/

ENV DEBIAN_FRONTEND=noninteractive
WORKDIR /opt

# -----------------------------
# 1) System dependencies
# -----------------------------
# - openjdk: needed for FastQC + Trimmomatic (both Java)
# - bowtie2, samtools, subread: alignment + BAM ops + featureCounts
# - python3-venv: for python venv if you ever want without uv (optional)
RUN apt-get update && apt-get install -y --no-install-recommends \
    ca-certificates \
    curl \
    unzip \
    gzip \
    pigz \
    perl \
    openjdk-17-jre-headless \
    python3 \
    python3-venv \
    bowtie2 \
    samtools \
    subread \
  && apt-get clean \
  && rm -rf /var/lib/apt/lists/*

# -----------------------------
# 2) FastQC install (Babraham)
# -----------------------------
RUN curl -LsSf https://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.12.1.zip \
      -o /tmp/fastqc.zip \
  && unzip /tmp/fastqc.zip -d /usr/local/src/ \
  && chmod +x /usr/local/src/FastQC/fastqc \
  && ln -sf /usr/local/src/FastQC/fastqc /usr/local/bin/fastqc \
  && rm -f /tmp/fastqc.zip

# -----------------------------
# 3) Trimmomatic install (Jar)
# -----------------------------
# Uses the official GitHub release jar.
# If you want a different version, change TRIMMOMATIC_VERSION.
ARG TRIMMOMATIC_VERSION=0.39
RUN mkdir -p /opt/trimmomatic \
  && curl -LsSf \
     https://github.com/usadellab/Trimmomatic/releases/download/v${TRIMMOMATIC_VERSION}/Trimmomatic-${TRIMMOMATIC_VERSION}.zip \
     -o /tmp/trimmomatic.zip \
  && unzip /tmp/trimmomatic.zip -d /opt/trimmomatic \
  && rm -f /tmp/trimmomatic.zip \
  && ln -sf /opt/trimmomatic/Trimmomatic-${TRIMMOMATIC_VERSION}/trimmomatic-${TRIMMOMATIC_VERSION}.jar /usr/local/bin/trimmomatic.jar

# Convenience wrapper so you can run `trimmomatic` directly.
RUN printf '%s\n' \
'#!/usr/bin/env bash' \
'set -euo pipefail' \
'java -jar /usr/local/bin/trimmomatic.jar "$@"' \
> /usr/local/bin/trimmomatic \
 && chmod +x /usr/local/bin/trimmomatic

# -----------------------------
# 4) Python tools: umi_tools + multiqc
# -----------------------------
# Create venv using uv and install python packages into it.
RUN uv venv /opt/venv
ENV VIRTUAL_ENV=/opt/venv
ENV PATH="/opt/venv/bin:/usr/local/bin:$PATH"

RUN uv pip install --no-cache-dir \
    multiqc \
    umi_tools

# -----------------------------
# 5) Verify installed tools
# -----------------------------
RUN for cmd in fastqc multiqc bowtie2 samtools featureCounts trimmomatic umi_tools; do \
      command -v "$cmd" >/dev/null 2>&1 || (echo "Missing $cmd" && exit 1); \
    done

# Default shell
ENTRYPOINT ["/bin/bash"]
