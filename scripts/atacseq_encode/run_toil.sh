#!/bin/sh
curdir=`pwd`
jobStore=$curdir/jobStore
cwltoil --restart \
    --no-container --jobStore $jobStore --workDir $curdir \
    --outdir $curdir \
    --beta-dependency-resolvers-configuration ~/work/micgent/python/lib/MICGENT/data/cwl_config/dep_resolver.yaml \
    --batchSystem Torque \
    --linkImports \
    ~/work/micgent/python/lib/MICGENT/data/cwl/atac_encode_samples.cwl \
    atac_encode_samples.yaml
