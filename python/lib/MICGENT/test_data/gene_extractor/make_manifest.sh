fastq_dir=$1
shift
[ -n "$fastq_dir ] | exit 1
python -m MICGENT.workflow_util iterate-sequence-run --from-dir $fastq_dir \
    --out-file-format csv --out-file sample_manifest.txt "[0-9]*.fastq.gz"
