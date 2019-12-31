#!/bin/sh
set -ex
cromwell -Dconfig.file=wdl.conf run ~/work/medi_scripts/16S/bbduk.wdl  bbduk_inputs.json bbduk_wf_opt.json
find wdl_final/ -name '*.fasta.gz' | grep 'collect_outputs/execution/glob'

