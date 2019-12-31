#!/bin/bash
set -ex
this_dir=$(cd $(dirname $0) && pwd)
curdir=$(pwd)

INPUT=$1; shift
[ -n "$INPUT" ] || exit 1

INPUT=$(cd "$(dirname "$INPUT")"; pwd)/$(basename "$INPUT")

[ -n "$CONDA_PREFIX_1" ] || exit 1
. $CONDA_PREFIX_1/etc/profile.d/conda.sh

conda activate atac-wdl

MICGENT_WDL=$this_dir
WDL=$MICGENT_WDL/atac_many.wdl
ATAC_HOME=$HOME/work/atac-seq-pipeline
ATAC_WDL=$ATAC_HOME/atac.wdl
BACKEND_CONF=$MICGENT_WDL/wdl.conf
WF_OPT=$MICGENT_WDL/atac_many_wf_opt.json
BACKEND=PBS

mkdir -p work && pushd work
export PATH="$ATAC_HOME/src:$PATH"

rm -f wdls.zip; zip -j wdls.zip $MICGENT_WDL/*.wdl $ATAC_WDL

cromwell -Dconfig.file=${BACKEND_CONF} -Dbackend.default=${BACKEND} run \
    ${WDL} \
    -i ${INPUT} \
    -o ${WF_OPT} \
    --imports wdls.zip \
    --metadata-output metadata.json \
    2>cromwell.err 1>cromwell.out
