#!/bin/bash
set -ex
this_dir=$(cd $(dirname $0) && pwd)
export PATH=$this_dir/..:$PATH
rm -f wdls.zip; zip -j wdls.zip $this_dir/*.wdl
cp $this_dir/../bbduk_adapters.fa ./
cromwell \
    -Dconfig.file=$this_dir/wdl.conf \
    run \
    $this_dir/dada.wdl  \
    $this_dir/dada_inputs.json \
    $this_dir/mtg16s_wf_opt.json \
    - \
    wdls.zip 2>cromwell.err 1>cromwell.out
#find wdl_final/ -name '*.fasta.gz' | grep 'collect_outputs/execution/glob'

