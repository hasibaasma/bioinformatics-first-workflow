#!/bin/bash
set -euo pipefail

# -----------------------------
# User settings
# -----------------------------
HOST_GENOME="genome.fa"        # FASTA used to build Bowtie2 index (if needed)
HOST_INDEX="host_reference"    # bowtie2 index prefix
THREADS=8

# Directories
LOGDIR="logs"
TRASHDIR="tmp"                 # temp files (keep or remove later)
MULTIQC_DIR="multiqc_report"

mkdir -p "${LOGDIR}"
mkdir -p "${TRASHDIR}"
mkdir -p "${MULTIQC_DIR}"

# -----------------------------
# 1. Build bowtie2 index (if not already done)
# -----------------------------
if [ ! -f "${HOST_INDEX}.1.bt2" ]; then
    echo "Building Bowtie2 index for ${HOST_GENOME}..." | tee "${LOGDIR}/bowtie2-build.log"
    bowtie2-build "${HOST_GENOME}" "${HOST_INDEX}" 2>&1 | tee -a "${LOGDIR}/bowtie2-build.log"
fi

# -----------------------------
# 2. Loop over samples with logging
# -----------------------------
for R1 in *_R1.fastq.gz; do
    SAMPLE=$(basename "$R1" _R1.fastq.gz)
    R2="${SAMPLE}_R2.fastq.gz"

    echo "==== Processing sample: ${SAMPLE} ====" | tee -a "${LOGDIR}/${SAMPLE}.master.log"

    # --- Step A: Host removal with bowtie2 + samtools ---
    # Capture bowtie2 stderr into bowtie2 log (mapping summary), capture samtools stderr to samtools log
    echo "[`date`] bowtie2 mapping start for ${SAMPLE}" | tee -a "${LOGDIR}/${SAMPLE}.master.log"

    # Run bowtie2, pipe to samtools; capture bowtie2 stderr to bowtie2 log, and samtools stderr to samtools log.
    # We still produce a BAM of UNMAPPED PAIRS (samtools view filtering), then sort by name.
    bowtie2 -x "${HOST_INDEX}" \
        -1 "${R1}" -2 "${R2}" \
        --very-sensitive-local -p "${THREADS}" \
        2> >(tee -a "${LOGDIR}/${SAMPLE}.bowtie2.log" >&2) \
        | samtools view -b -f 12 -F 256 - \
          2> >(tee -a "${LOGDIR}/${SAMPLE}.samtools_view.log" >&2) \
        | samtools sort -n -o "${TRASHDIR}/${SAMPLE}_host_removed.bam" \
          2> >(tee -a "${LOGDIR}/${SAMPLE}.samtools_sort.log" >&2)

    echo "[`date`] bowtie2 mapping finished for ${SAMPLE}" | tee -a "${LOGDIR}/${SAMPLE}.master.log"

    # Flagstat on filtered BAM (log)
    samtools flagstat "${TRASHDIR}/${SAMPLE}_host_removed.bam" \
        2> >(tee -a "${LOGDIR}/${SAMPLE}.samtools_flagstat.log" >&2) \
        | tee -a "${LOGDIR}/${SAMPLE}.samtools_flagstat.log"

    # --- Step B: Convert BAM -> paired FASTQ (host-removed) ---
    echo "[`date`] samtools fastq (BAM->FASTQ) for ${SAMPLE}" | tee -a "${LOGDIR}/${SAMPLE}.master.log"

    samtools fastq \
        -1 "${TRASHDIR}/${SAMPLE}_hostRemoved_R1.fastq.gz" \
        -2 "${TRASHDIR}/${SAMPLE}_hostRemoved_R2.fastq.gz" \
        -0 /dev/null -s /dev/null -n \
        "${TRASHDIR}/${SAMPLE}_host_removed.bam" \
        2> >(tee -a "${LOGDIR}/${SAMPLE}.samtools_fastq.log" >&2)

    # --- Step C: Quality trimming with fastp ---
    echo "[`date`] fastp trimming for ${SAMPLE}" | tee -a "${LOGDIR}/${SAMPLE}.master.log"

    fastp \
      --in1 "${TRASHDIR}/${SAMPLE}_hostRemoved_R1.fastq.gz" \
      --in2 "${TRASHDIR}/${SAMPLE}_hostRemoved_R2.fastq.gz" \
      --out1 "${SAMPLE}_trimmed_R1.fastq.gz" \
      --out2 "${SAMPLE}_trimmed_R2.fastq.gz" \
      --cut_right --cut_window_size 4 --cut_mean_quality 20 \
      -l 50 \
      --detect_adapter_for_pe \
      -y \
      --thread "${THREADS}" \
      --html "${SAMPLE}_fastp.html" \
      --json "${SAMPLE}_fastp.json" \
      2>&1 | tee -a "${LOGDIR}/${SAMPLE}.fastp.log"

    echo "[`date`] fastp finished for ${SAMPLE}" | tee -a "${LOGDIR}/${SAMPLE}.master.log"

    # Optionally remove intermediate host-removed FASTQs / BAM (uncomment to save space)
    # rm -f "${TRASHDIR}/${SAMPLE}_host_removed.bam" "${TRASHDIR}/${SAMPLE}_hostRemoved_R1.fastq.gz" "${TRASHDIR}/${SAMPLE}_hostRemoved_R2.fastq.gz"

    echo "==== Finished sample: ${SAMPLE} ====" | tee -a "${LOGDIR}/${SAMPLE}.master.log"
done

# -----------------------------
# 3. MultiQC report (capture output)
# -----------------------------
echo "[`date`] Running MultiQC ..." | tee -a "${LOGDIR}/multiqc.master.log"
multiqc . -o "${MULTIQC_DIR}" 2>&1 | tee -a "${LOGDIR}/multiqc.master.log"
echo "[`date`] MultiQC finished" | tee -a "${LOGDIR}/multiqc.master.log"

echo "✅ All done. Logs are in the '${LOGDIR}' folder. MultiQC report: ${MULTIQC_DIR}/multiqc_report.html"
