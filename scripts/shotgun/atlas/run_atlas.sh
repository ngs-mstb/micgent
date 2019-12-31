#!/bin/sh
atlas assemble --jobs 500 --out-dir final config.yaml \
    --cluster '"export SHELL=/bin/sh; qsub -V -S /bin/sh -l nodes=1:ppn={threads}"' \
    --latency-wait 60 \
    --keep-going \
    --printshellcmds --verbose 

