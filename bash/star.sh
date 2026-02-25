# ls -la /work/natanastas/Merge/20250127 | awk '{print $9}' | awk -F "_L" '{print $1}' | sort | uniq > SampleList

pathToFASTQFiles="$1"
SampleList="$2"

pathToGenomeRef="$3"
pathToGenomeAnnotation="$4"

numberOfThreads="$5"


printf "mkdir quality\n"
printf "mkdir quality/raw\n"
printf "mkdir quality/trim\n"
printf "\n\n"

printf "source ~/.bashrc\n"
printf "\n\n"

printf "conda activate fastqc_v0.12.1"
printf "\n\n"

cat $SampleList | while read line; do

	printf "### $line ###\n"

	printf "fastqc \
	-t $numberOfThreads \
	$pathToFASTQFiles/$line*R1*.fastq.gz \
 	$pathToFASTQFiles/$line*R2*.fastq.gz \
 	-o ./quality/raw/ \
 	\n"

 	printf "\n"
 done;

 printf "conda deactivate"
 printf "\n\n"

 printf "conda activate trim_galore"
 printf "\n\n"

 cat $SampleList | while read line; do

	printf "### $line ###\n"

	printf "trim_galore \
	--path_to_cutadapt /apps/miniconda3/envs/trim_galore/bin/cutadapt \
	--length 50 \
	--quality 28 \
	--fastqc \
	--paired \
	$pathToFASTQFiles/$line*R1*.fastq.gz \
	$pathToFASTQFiles/$line*R2*.fastq.gz \
	-o ./quality/trim/ \
	\n"

	printf "\n"
 done;

printf "conda deactivate"
printf "\n\n"

printf "conda activate star_aligner"
printf "\n\n"

cat $SampleList | while read line; do

	printf "### $line ###\n"

	printf "STAR \
	--genomeDir $pathToGenomeRef \
	--runThreadN $numberOfThreads \
	--readFilesCommand zcat \
	--readFilesIn \
	./quality/trim/$line*R1*val_1.fq.gz \
	./quality/trim/$line*R2*val_2.fq.gz \
	--outFileNamePrefix $line. \
	--outReadsUnmapped Fastx \
    \n"

	printf "samtools view     -@ $numberOfThreads -bS $line.Aligned.out.sam -o $line.bam \n"
	printf "rm $line.Aligned.out.sam\n"
	printf "\n"
done;

printf "conda deactivate"
printf "\n\n"


cat $SampleList | while read line; do

	printf "### $line ###\n"

	printf "samtools flagstat -@ $numberOfThreads $line.bam > $line.report.txt \n"
	printf "samtools view 	  -@ $numberOfThreads -b -F 4 $line.bam -o $line.mapped.bam \n"
	printf "samtools sort     -@ $numberOfThreads $line.mapped.bam -o $line.mapped.sorted.bam \n"

	printf "\n"
done;

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
