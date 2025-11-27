
pathToFASTQFiles="$1"
SampleList="$2"

cat $SampleList | while read line; do

	printf "### $line ###\n"

	printf "zcat "
	printf $pathToFASTQFiles/$line
	printf "*R1*.fastq.gz > "
	printf $line
	printf "_L001_R1_001.fastq\n"

    printf "zcat "
	printf $pathToFASTQFiles/$line
	printf "*R2*.fastq.gz > "
	printf $line
	printf "_L001_R2_001.fastq\n"

	printf "gzip "
	printf $line
	printf "_L001_R1_001.fastq\n"

	printf "gzip "
	printf $line
	printf "_L001_R2_001.fastq\n"

	printf "\n\n"
done;
