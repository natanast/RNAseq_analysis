# How to Use the RNA-seq Scripts

## STAR Alignment Script
This script runs a full RNA-seq preprocessing and alignment pipeline using:

- FastQC (for quality control),

- Trim Galore (for trimming adapters and low-quality bases),

- STAR (for genome alignment),

- SAMtools (for post-alignment processing), and

- featureCounts (for generating gene-level read counts).

##  Run the Pipeline

```bash
bash runSTAR.sh \
/path_to_FASTQ_files \
SampleList \
/work_1/nikospech/hg38/gencode/star \
/work_1/nikospech/hg38/gencode/gencode.v47.primary_assembly.basic.annotation.gtf \
16
```

## Output

./quality/raw/: Raw FastQC reports

./quality/trim/: Trimmed FASTQ files + FastQC post-trim reports

*.bam, *.mapped.bam, *.mapped.sorted.bam: Intermediate BAM files

gene-counts.txt: Final gene-level count matrix

*.report.txt: Alignment statistics