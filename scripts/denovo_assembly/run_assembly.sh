#!/bin/sh
export PATH="/Xxx/ID/other_install/miniconda2/bin:/Xxx/ID/other_install/bin:$PATH"
source activate assembly
set -ex
MICGENT_DATA=/Xxx/xxx/micgent_db python -m MICGENT.denovo_mic_assembly assemble-de-novo
./assembly_batch.mkf.bat
python -m MICGENT.denovo_mic_assembly collect-bugbuilder-asm-metrics "assembly/*/out.log" asm_metrics.csv
python -m MICGENT.denovo_mic_assembly export-assembly "assembly/*/out" asm_export

