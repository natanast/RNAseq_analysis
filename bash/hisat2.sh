pathToFASTQFiles="/work_1/nikospech/20240611_mergeLanes_RNAseq_Glykeria/clean-data/"

printf "mkdir qualityRaw\n"

cat SampleList | while read line; do

	printf "### Sample: "
	printf $line
	printf "  ### \n"

	printf "fastqc "
	printf $pathToFASTQFiles
	printf $line
	printf "_R1_001.fastq.gz "
	printf $pathToFASTQFiles
	printf $line
	printf "_R2_001.fastq.gz "
	printf " -o qualityRaw/"
	printf "\n"
	
	printf "/home/bio_tmp/HumanNGS/Apps/trim_galore/trim_galore"
	printf " --path_to_cutadapt /home/bio_tmp/HumanNGS/Apps/cutadapt-1.8.1/bin/cutadapt"
	printf " --length 50 --quality 28 --fastqc --paired "
	printf $pathToFASTQFiles
	printf $line
	printf "_R1_001.fastq.gz "
	printf $pathToFASTQFiles
	printf $line
	printf "_R2_001.fastq.gz"
	printf " -o ."
	printf "\n"

	printf "hisat2 --threads 16"
	printf " -x /work_1/nikospech/hg38/hisat2/genome"
	printf " -1 ./"
	printf $line
	printf "_R1_001_val_1.fq.gz "
	printf " -2 ./"
	printf $line
	printf "_R2_001_val_2.fq.gz "
	printf " -S "
	printf $line
	printf ".sam"
	printf " --un-conc-gz "
	printf $line
	printf ".unmap.fastq.gz"
	printf "\n"

	printf "fastqc "
	printf $line
	printf ".unmap.fastq.1.gz "
	printf $line
	printf ".unmap.fastq.2.gz"
	printf "\n"

	printf "/home/bio_tmp/HumanNGS/Apps/trim_galore/trim_galore"
	printf " --path_to_cutadapt /home/bio_tmp/HumanNGS/Apps/cutadapt-1.8.1/bin/cutadapt"
	printf " --length 50 --three_prime_clip_R1 75 --three_prime_clip_R2 75 --quality 28 --fastqc --paired "
	printf $line
	printf ".unmap.fastq.1.gz "
	printf $line
	printf ".unmap.fastq.2.gz"
	printf " -o ."
	printf "\n"

	printf "hisat2 --threads 16"
	printf " -x /work_1/nikospech/hg38/hisat2/genome"
	printf " -1 ./"
	printf $line
	printf ".unmap.fastq.1.gz_val_1.fq.gz"
	printf " -2 ./"
	printf $line
	printf ".unmap.fastq.2.gz_val_2.fq.gz"
	printf " -S "
	printf $line
	printf ".p2.sam"
	printf " --un-conc-gz "
	printf $line
	printf ".unmap2.fastq.gz"
	printf "\n"

	printf "fastqc "
	printf $line
	printf ".unmap2.fastq.1.gz "
    printf $line
	printf ".unmap2.fastq.2.gz"
    printf "\n"

	printf "samtools view -@ 16 -bS "
	printf $line
	printf ".sam -o "
	printf $line
	printf ".bam"
	printf "\n"

	printf "samtools view -@ 16 -bS "
	printf $line
	printf ".p2.sam -o "
	printf $line
	printf ".p2.bam"
	printf "\n"

	printf "samtools flagstat -@ 16 "
	printf $line
	printf ".bam > "
	printf $line
	printf ".report.txt"
	printf "\n"

	printf "samtools flagstat -@ 16 "
	printf $line
	printf ".p2.bam > "
	printf $line
	printf ".p2.report.txt"
	printf "\n"

	printf "samtools view -@ 16 -b -F 4 "
	printf $line
	printf ".bam -o "
	printf $line
	printf ".mapped.bam"
	printf "\n"

	printf "samtools view -@ 16 -b -F 4 "
	printf $line
	printf ".p2.bam -o "
	printf $line
	printf ".p2.mapped.bam"
	printf "\n"

	printf "samtools sort -@ 16 "
	printf $line
	printf ".mapped.bam -o "
	printf $line
	printf ".mapped.sorted.bam"
	printf "\n"

	printf "samtools sort -@ 16 "
	printf $line
	printf ".p2.mapped.bam -o "
	printf $line
	printf ".p2.mapped.sorted.bam"
	printf "\n"

	printf "samtools merge -@ 16 "
	printf $line
	printf ".merged.bam "
	printf $line
	printf ".mapped.sorted.bam "
	printf $line
	printf ".p2.mapped.sorted.bam "
	printf "\n"

	printf "samtools sort -@ 16 "
	printf $line
	printf ".merged.bam -o "
	printf $line
	printf ".merged.sorted.bam"
	printf "\n"

	printf "\n\n"

	done;
	
printf "rm *.sam"
printf "\n\n"

printf "/home/bio_tmp/HumanNGS/Apps/subread-1.6.4-Linux-x86_64/bin/featureCounts"
printf " -a /work_1/nikospech/hg38/genomic.gff"
printf " -T 16 "
printf " -t gene "
printf " -g ID "
printf " -p"
printf " -o gene-counts.txt "
printf "*.merged.sorted.bam"
printf "\n"
