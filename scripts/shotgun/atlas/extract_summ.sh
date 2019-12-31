#!/bin/sh
python -m MICGENT.atlas extract-func-summary-many \
    --samp-id-extractor '^[^_]+' \
    "SP*/megahit_21_121_20_normalization_k21_t100/SP*_annotations.txt" func.txt
