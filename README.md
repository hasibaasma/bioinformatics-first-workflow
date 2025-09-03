

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

### Step 1 – Quality Control (example)

**Purpose:** Remove low-quality reads and adapter sequences
**Tools:**   `Bowtie2`, `samtools`
**Inputs:** Subsampled FASTQ files (from `data/fastq_subsampled/`)
**Outputs:** Cleaned FASTQs, QC reports (`.html`, `.json`, or `.txt`)
**Command:**

```bash
fastp --in1 sample.fastq.gz --out1 cleaned.fastq.gz ...
```

---

### Step 2 ...

**Purpose:** ...
**Tools:** ...
**Inputs:** ...
**Outputs:** ...
**Command:**


---

### Step X – Analysis (e.g., DESeq2, variant calling, etc.)

**Purpose:** ...
**Tools:** ...
**Inputs:** ...
**Outputs:** ...
**Command:**

