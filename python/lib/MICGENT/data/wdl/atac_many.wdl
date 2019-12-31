import "atac.wdl" as atac_wf

## Run ENCODE WDL ATAC-Seq pipeline on multiple samples
## Output file types can be found here (docs for bigscript version):
## https://github.com/kundajelab/atac_dnase_pipelines
workflow atac_many {
  Array[Pair[String,Array[Array[Array[File]]]]] sample_fastqs # samp_id => [rep_id][merge_id][read_end_id]
  scatter(sample in sample_fastqs) {
    String sample_id = sample.left
    call atac_wf.atac as atac {
      input:
      fastqs = sample.right
    }
    call label_meta as label_sample_id {
      input:
      label = sample_id
    }
    call collect_outputs as col_bams {
      input:
      inp_files = atac.out_nodup_bams,
      prefix = sample_id
    }
    call collect_outputs as col_tas {
      input:
      inp_files = atac.out_tas,
      prefix = sample_id
    }
    call collect_outputs as col_bfilt_npeaks {
      input:
      inp_files = atac.out_bfilt_npeaks,
      prefix = sample_id
    }
    call collect_outputs as col_qc_report {
      input:
      inp_files = atac.out_qc_report,
      prefix = sample_id
    }
    call collect_outputs as col_qc_json {
      input:
      inp_files = atac.out_qc_json,
      prefix = sample_id
    }
  }
  output {
    Array[Pair[String,Array[File]]] out_bams = col_bams.out_files
    Array[Pair[String,Array[File]]] out_tas = col_tas.out_files
    Array[Pair[String,Array[File]]] out_bfilt_npeaks = col_bfilt_npeaks.out_files
    Array[Pair[String,Array[File]]] out_qc_report = col_qc_report.out_files
    Array[Pair[String,Array[File]]] out_qc_json = col_qc_json.out_files
  }
}

## A noop task to leave a labeled node
## in Cromwell output metadata file (cromwell --metadata-output).
## Such labeled nodes can help querying metadata JSON such as
## extracting specific deeply nested task outputs from scatter
## iterations (e.g. use http://jmespath.org/ projection queries
## to find labeled node and then go to its container and get sister
## node with scatter iteration output)
task label_meta {
  String label
  command {}
  output {
    String out_label = label
  }
}

task collect_outputs {
  Array[File] inp_files
  String prefix=""
  String prefix_sep="_"
  String prefix_fn = if prefix!="" then prefix + prefix_sep else prefix
  command {
    mkdir final
    cat ${write_lines(inp_files)} | \
    xargs -I @@ sh -c 'cp -l @@ final/${prefix_fn}$(basename @@)'
  }
  output {
    Pair[String,Array[File]] out_files=(prefix,glob("final/*"))
  }
}

