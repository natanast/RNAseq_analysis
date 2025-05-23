pathToFASTQFiles="$1"
SampleList="$2"

pathToGenomeRef="$3"
pathToGenomeAnnotation="$4"

numberOfThreads="$5"

bash star.sh $pathToFASTQFiles $SampleList $pathToGenomeRef $pathToGenomeAnnotation $numberOfThreads > run.sh
bash run.sh