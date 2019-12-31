import "bbduk.wdl" as bbduk_wf
import "subsample_paired_reads.wdl" as subsamp_wf

workflow dada_wf {
  File sampleSheetFile
  Array[String] primer_literals
  File? adapter_file
  Int threads
  Float? samplerate
  Int? max_reads
  Array[Object] inputSamples = read_objects(sampleSheetFile)
  call bbduk_wf.bbduk_cut_16s_primers as cleanup {
    input: inputSamples=inputSamples,
    primer_literals=primer_literals,
    adapter_file=adapter_file,
    bbduk_threads=threads,
    samplerate=samplerate,
    reads=max_reads
    }
  #call subsamp_wf.subsample_paired_reads_collection as subsamp {
  #  input: SampleID = cleanup.out_SampleID,
  #  inp_seqs = cleanup.out_seqs
  #}
  call dada {
    input: SampleID = cleanup.out_SampleID,
    inp_seqs = cleanup.out_seqs,
    threads = threads 
  }
  output {
    File tab = dada.tab
    File otu = dada.otu
    File taxonomy = dada.taxonomy
    File report = dada.report
  }
}

task dada {
    Array[String] SampleID
    Array[Array[File]] inp_seqs
    Int threads=4
    Array[String] inp_dd_SampleID = SampleID

    command <<<
      set -ex
      python -m MICGENT.workflow_util join-files-with-ids \
      ${write_tsv(inp_dd_SampleID)} \
      ${write_tsv(inp_seqs)} \
      sample_sheet.tsv
      dada_pipeline.R \
        --sample_sheet_file sample_sheet.tsv \
        --report_seq_quality \
        --threads ${threads}
      tar -czf dada_report.tgz *.html plots data widget_deps 
    >>>
    output {
      File tab = "tab.csv"
      File otu = "otu.shared"
      File taxonomy = "cons.taxonomy"
      File report = "dada_report.tgz"
    }
    runtime {
        cpu : threads
    }
}

task collect_outputs {
  Array[Array[File]] inp_seqs
  command {
    mkdir final
    cat ${write_lines(inp_seqs)} | tr '\t' '\n' | \
    xargs -I @@ cp @@ final/
  }
  output {
    Array[File] out_seqs=glob("final/*")
  }
}

