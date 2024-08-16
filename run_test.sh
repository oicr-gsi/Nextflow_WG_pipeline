module load biomodal

ts=`date '+%Y-%m-%d.%H.%M.%S'`
log="logs/nextflow.$ts.log"


nextflow -log $log run \
  ./main.nf \
  -resume 
