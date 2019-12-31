#!/bin/sh
cwltool --debug --disable-pull --no-container --cachedir `pwd` \
    --beta-dependency-resolvers-configuration ~/work/micgent/python/lib/MICGENT/data/cwl_config/dep_resolver.yaml \
    ~/work/micgent/python/lib/MICGENT/data/cwl/atac_encode_samples.cwl \
    atac_encode_samples.yaml

