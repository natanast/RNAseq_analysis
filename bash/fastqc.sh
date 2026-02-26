# ls -la /work/natanastas/Merge/20250127 | awk '{print $9}' | awk -F "_L" '{print $1}' | sort | uniq > SampleList

pathToFASTQFiles="$1"
SampleList="$2"

numberOfThreads="$3"

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