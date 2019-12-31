workflow bbduk_cut_16s_primers {
  Array[Object] inputSamples
  Array[String] primer_literals
  File? adapter_file
  Int bbduk_threads
  Float? samplerate
  Int? reads
  scatter (sample in inputSamples) {
  call bbduk_cut_primers as cut_left {
    input: SampleID=sample.SampleID,
    inp_seqs=[sample.file1,sample.file2],
    adapter_literals=primer_literals,
    ktrim='l',
    rcomp='f',
    threads=bbduk_threads,
    samplerate=samplerate,
    reads=reads
    }
  call bbduk_cut_primers as cut_right {
    input: SampleID=cut_left.out_SampleID,
    inp_seqs=cut_left.out_seqs,
    adapter_literals=primer_literals,
    adapter_file=adapter_file,
    ktrim='r',
    rcomp='t',
    threads=bbduk_threads
    }
  }
  output {
    Array[String] out_SampleID = cut_right.out_SampleID
    Array[Array[File]] out_seqs=cut_right.out_seqs
    Array[Array[File]] out_stats=transpose([cut_left.out_stats,cut_right.out_stats])    
  }
}

task bbduk_cut_primers {
    String SampleID
    Int threads=4
    Array[File] inp_seqs
    Array[String]? adapter_literals
    File? adapter_file
    String ktrim
    String rcomp
    Float? samplerate
    Int? reads
    String out_stats_fn="${SampleID}.bbduk.stats"

    command <<<
        adapter_str="${sep=',' adapter_literals}"
        [ -n "$adapter_str" ] && adapter_str="literal=$adapter_str"
        bbduk.sh \
        in=${inp_seqs[0]} in2=${inp_seqs[1]} \
        out=${SampleID}_1.fastq.gz out2=${SampleID}_2.fastq.gz  \
        $adapter_str \
        ${'ref='+adapter_file} \
        ktrim=${ktrim} k=13 hdist=1 \
        rcomp=${rcomp} \
        ${"samplerate="+samplerate} \
        ${"reads="+reads} \
        tbo tpe \
        stats=${out_stats_fn} \
        ${"threads="+threads} \
        copyundefined
    >>>
    output {
        Array[File] out_seqs = ["${SampleID}_1.fastq.gz", 
            "${SampleID}_2.fastq.gz"]
        String out_SampleID="${SampleID}"
        File out_stats="${out_stats_fn}"
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

