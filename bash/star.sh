# ls -la /work/natanastas/Merge/20250127 | awk '{print $9}' | awk -F "_L" '{print $1}' | sort | uniq > SampleList

pathToGenomeRef="/work_1/nikospech/hg38/gencode/star"
pathToGenomeAnnotation="/work_1/nikospech/hg38/gencode/gencode.v47.primary_assembly.basic.annotation.gtf"
pathToFASTQFiles="/mnt/new_home/kate_mallou/Karolinska_RNASeq/"
numberOfThreads=16

printf "mkdir quality"
printf "mkdir quality/raw\n"
printf "mkdir quality/trim\n"

printf "source activate star_aligner"
printf "\n\n"

cat SampleList | while read line; do

	printf "### Sample: "
	printf $line
	printf "  ### \n"

	printf "fastqc -t "
    printf $numberOfThreads
    printf " "
	printf $pathToFASTQFiles
	printf $line
	printf "*_R1_*.fastq.gz "
	printf $pathToFASTQFiles
	printf $line
	printf "*_R2_*.fastq.gz "
	printf " -o ./quality/raw/"
	printf "\n"

	printf "/home/bio_tmp/HumanNGS/Apps/trim_galore/trim_galore"
	printf " --path_to_cutadapt /home/bio_tmp/HumanNGS/Apps/cutadapt-1.8.1/bin/cutadapt"
	printf " --length 50 --quality 28 --fastqc --paired  "
	printf $pathToFASTQFiles
	printf $line
	printf "*_R1_*.fastq.gz "
	printf $pathToFASTQFiles
	printf $line
	printf "*_R2_*.fastq.gz"
	printf " -o ./quality/trim/"
	printf "\n"

	printf "STAR"
	printf " --genomeDir "
    printf $pathToGenomeRef
   	printf " --runThreadN "
    printf $numberOfThreads
    printf " --readFilesCommand zcat"
    printf " --readFilesIn ./quality/trim/"
    printf $line
    printf "*_R1_*val_1.fq.gz ./quality/trim/"
    printf $line
    printf "*_R2_*val_2.fq.gz"
    printf " --outFileNamePrefix "
    printf $line
	printf "."
	printf " --outReadsUnmapped Fastx"
    printf "\n"

	printf "samtools view"
	printf " -@ "
	printf $numberOfThreads
	printf " -bS "
	printf $line
	printf ".Aligned.out.sam -o "
	printf $line
	printf ".bam"
	printf "\n"

	printf "samtools flagstat"
	printf " -@ "
	printf $numberOfThreads
	printf " "
	printf $line
	printf ".bam > "
	printf $line
	printf ".report.txt"
	printf "\n"

	printf "samtools view"
	printf " -@ "
	printf $numberOfThreads
	printf " -b -F 4 "
	printf $line
	printf ".bam -o "
	printf $line
	printf ".mapped.bam"
	printf "\n"

	printf "samtools sort"
	printf " -@ "
	printf $numberOfThreads
	printf " "
	printf $line
	printf ".mapped.bam -o "
	printf $line
	printf ".mapped.sorted.bam"
	printf "\n"

	printf "\n\n"
done;

printf "conda deactivate"
printf "\n\n"
printf "rm *.sam"
printf "\n\n"

printf "/home/bio_tmp/HumanNGS/Apps/subread-1.6.4-Linux-x86_64/bin/featureCounts"
printf " -a "
printf $pathToGenomeAnnotation
printf " -T "
printf $numberOfThreads
printf " -t gene "
printf " -g gene_name "
printf " -p"
printf " -o gene-counts.txt "
printf "*.mapped.sorted.bam"
printf "\n"
