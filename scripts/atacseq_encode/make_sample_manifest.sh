#!/bin/sh
set -ex
python -m MICGENT.workflow_util iterate-sequence-run \
    --samp-id-extractor '^xxx-([^_]+)_' \
    --sample-id-prefix AT \
    --out-file samples.csv --out-file-format csv \
    --allow-repeated-sample-id "FASTQ/*/xxx-*.fastq.gz"
python -m MICGENT.workflow_util manifest-to-atacseq-cwl \
    samples.csv samples.yaml \
    '^(.*)[ab]'
