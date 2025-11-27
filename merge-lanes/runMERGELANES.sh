pathToFASTQFiles="$1"
SampleList="$2"

bash mergeLanes.sh $pathToFASTQFiles $SampleList > run.sh
bash run.sh