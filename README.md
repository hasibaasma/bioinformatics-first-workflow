

#  README

## 1. Data sources

This task is based on publicly available data from a study of **\Dynamics of Gut Microbiota After Fecal Microbiota Transplantation in Ulcerative Colitis: Success Linked to Control of Prevotellaceae**. The dataset includes multiple samples under different conditions (e.g., treated vs. control) and was originally sequenced using **\Illumina NovaSeq platform (100 bp single end read)**.
The subsampled and cleaned FASTQs are stored in `data/` and are used as the inputs for the workflow.

---

## 2. How to download

INSTRUCTIONS TO ACCESS THE DATA
### Example using SRA Toolkit

```bash
CODE TO DOWNLOAD
```


---

## 3. Pre-processing / subsampling

INCLUDE THE METHOD YOU USED TO SUBSAMPLE, MINATURIZE, OR TRIM DOWN

1. **STEP 1** ...

Example:

```bash
CODE TO SUBSAMPLE
```


---

## 4. How the workflow works
Removal of host genome uses **\Bowtie2 to map reads to reference genome
Filter unmapped paired reads from output using **\Samtools
Quality Trimming, adapter removal and low complexity filtering is done with **\ Fastp
Reports of mapping and quality trimming are generated with **\ MultiQC
The workflow files is stored in `workflow/`.

---

### Step 1 – Host Removal

**Purpose:** Remove host-derived reads (e.g., human DNA) to focus on microbial content
**Tools:**   Bowtie2, samtools
**Inputs:** Raw paired-end FASTQ files of subsample (*_R1.fastq.gz, *_R2.fastq.gz), host genome index
**Outputs:** Cleaned FASTQs, QC reports (`.html`, `.json`, or `.txt`)
**Command:**

```bash
fastp --in1 sample.fastq.gz --out1 cleaned.fastq.gz ...
```

---

### Step 2 Quality Trimming & Adapter Removal

**Purpose:** Remove low-quality bases and adapter sequences to improve downstream analyses
**Tools:** fastp
**Inputs:** Host-filtered paired-end FASTQ files (*_hostRemoved_R1.fastq.gz, *_hostRemoved_R2.fastq.gz)
**Outputs:** Cleaned paired-end FASTQ files (*_trimmed_R1.fastq.gz, *_trimmed_R2.fastq.gz), QC reports (*_fastp.html, *_fastp.json)
**Command:**

```bash
fastp \
  --in1 sample_hostRemoved_R1.fastq.gz \
  --in2 sample_hostRemoved_R2.fastq.gz \
  --out1 sample_trimmed_R1.fastq.gz \
  --out2 sample_trimmed_R2.fastq.gz \
  --cut_right --cut_window_size 4 --cut_mean_quality 20 \
  -l 50 \
  --detect_adapter_for_pe \
  -y \
  --thread 8 \
  --html sample_fastp.html \
  --json sample_fastp.json
```

---

### Step 3 – Generate Summary Report

**Purpose:** Aggregate QC results across all samples for easy interpretation
**Tools:** MultiQC
**Inputs:** Logs and reports from bowtie2, samtools, and fastp (*.html, *.json, *.txt)
**Outputs:** Combined QC report (multiqc_report.html)
**Command:**

```
multiqc . -o multiqc_report
```
