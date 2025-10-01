
pathToBamFiles="data/"
pathToCoordinates="bcrtcr.fa"
pathToIMGTReference="IMGT+C.fa"

numberOfThreads=4

printf "mkdir trust4\n"
printf "\n\n"

# printf "source activate trust4"
# printf "\n\n"

cat $SampleList | while read line; do

	printf "### $line ###\n"
	
	printf "run-trust4 \
	-t $numberOfThreads \
	-b $pathToBamFiles/$line.mapped.sorted.bam \
	-f $pathToCoordinates \
	--ref $pathToIMGTReference \
	--od ./trust4/$line \
	\n"
	
	printf "gzip ./trust4/$line/*.fa \n"
	printf "gzip ./trust4/$line/*.fq \n"

	printf "\n\n"
done;

printf "conda deactivate"