pathToFASTQFiles="$1"
SampleList="$2"

numberOfThreads="$3"

bash fastqc.sh $pathToFASTQFiles $SampleList $numberOfThreads > run.sh
bash run.sh