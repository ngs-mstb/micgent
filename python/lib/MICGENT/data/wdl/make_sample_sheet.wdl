task make_sample_sheet {
  String read_files_glob
  String samp_id_extractor
  String samp_sheet="samples.csv"
  Array[File] read_files
  command {
    python -m MICGENT.workflow_util iterate-sequence-run \
    --samp-id-extractor ${samp_id_extractor} \
    --out-file ${samp_sheet} \
    --out-file-format csv \
    "${read_files_glob}"
  }
  output {
    File samp_sheet=${samp_sheet}
  }
  runtime {
  conda_env: assembly
  }
}

