#/bin/sh
echo '(hostname; date; pwd; (printenv | sort)) | tee test.qsub.local.out' > test.qsub && qsub -S /bin/sh -d `pwd` -e `pwd`/test.qsub.err -o `pwd`/test.qsub.out test.qsub

