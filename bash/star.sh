# ls -la /work/natanastas/Merge/20250127 | awk '{print $9}' | awk -F "_L" '{print $1}' | sort | uniq > SampleList

pathToGenomeRef="/work_1/nikospech/hg38/gencode/star"
pathToGenomeAnnotation="/work_1/nikospech/hg38/gencode/gencode.v47.primary_assembly.basic.annotation.gtf"
pathToFASTQFiles="/mnt/new_home/kate_mallou/Karolinska_RNASeq/"
numberOfThreads=16

printf "mkdir quality\n"
printf "mkdir quality/raw\n"
printf "mkdir quality/trim\n"
printf "\n\n"

printf "source activate star_aligner"
printf "\n\n"

cat SampleList | while read line; do

	printf "### $line ###\n"

	printf "fastqc \
	-t $numberOfThreads \
	$pathToFASTQFiles/$line*_R1_*.fastq.gz \
	$pathToFASTQFiles/$line*_R2_*.fastq.gz \
	-o ./quality/raw/ \
	\n"

	printf "/home/bio_tmp/HumanNGS/Apps/trim_galore/trim_galore \
	--path_to_cutadapt /home/bio_tmp/HumanNGS/Apps/cutadapt-1.8.1/bin/cutadapt \
	--length 50 \
	--quality 28 \
	--fastqc \
	--paired \
	$pathToFASTQFiles/$line*_R1_*.fastq.gz \
	$pathToFASTQFiles/$line*_R2_*.fastq.gz \
	-o ./quality/trim/ \
	\n"

	printf "STAR \
	--genomeDir $pathToGenomeRef \
	--runThreadN $numberOfThreads \
	--readFilesCommand zcat \
	--readFilesIn \
	./quality/trim/$line*_R1_*val_1.fq.gz \
	./quality/trim/$line*_R2_*val_2.fq.gz \
	--outFileNamePrefix $line. \
	--outReadsUnmapped Fastx \
    \n"

	printf "samtools view \
	-@ $numberOfThreads \
	-bS $line.Aligned.out.sam \
	-o $line.bam \
	\n"

	printf "samtools flagstat \
	-@ $numberOfThreads \
	$line.bam > $line.report.txt \
	\n"

	printf "samtools view \
	-@ $numberOfThreads \
	-b -F 4 \
	$line.bam \
	-o $line.mapped.bam \
	\n"

	printf "samtools sort \
	-@ $numberOfThreads \
	$line.mapped.bam \
	-o $line.mapped.sorted.bam \
	\n"

	printf "\n\n"
done;

printf "conda deactivate"
printf "\n\n"

printf "rm *.sam"
printf "\n\n"

printf "/home/bio_tmp/HumanNGS/Apps/subread-1.6.4-Linux-x86_64/bin/featureCounts \
-a $pathToGenomeAnnotation \
-T $numberOfThreads \
-t gene \
-g gene_name \
-p \
-o gene-counts.txt \
*.mapped.sorted.bam \
\n"
