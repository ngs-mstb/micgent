#!/bin/bash
wf_cwl=$1
shift
wf_inp=$1
[ -n "$wf_cwl" ] | exit 1
[ -n "$wf_inp" ] | exit 1
. /xxxscr/xxx/mc3/etc/profile.d/conda.sh
conda activate toil
cwltool --leave-outputs --leave-tmpdir --debug \
    --beta-dependency-resolvers-configuration dep_resolver.yaml \
    $wf_cwl $wf_inp > cwltool.log 2>&1

