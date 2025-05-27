# How to Use the RNA-seq Scripts

## STAR Alignment Script
This command runs the scripts which uses STAR alignment for a batch of paired-end RNA-seq samples.

```bash
bash runSTAR.sh \
/path_to_FASTQ_files /
SampleList /
/work_1/nikospech/hg38/gencode/star /
/work_1/nikospech/hg38/gencode/gencode.v47.primary_assembly.basic.annotation.gtf /
16
```
