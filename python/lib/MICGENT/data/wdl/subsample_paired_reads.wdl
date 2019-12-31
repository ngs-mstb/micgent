
workflow subsample_paired_reads_collection {
  Array[String] SampleID
  Array[Array[File]] inp_seqs
  Float read_count
  String seed=1
  Array[String] inp_sprc_SampleID=SampleID
  Array[Pair[String,Array[File]]] inp_recs = zip(inp_sprc_SampleID,inp_seqs)
  
  scatter(inp_seq_rec in inp_recs) {
    call subsample_paired_reads as spr {
      input: SampleID=inp_seq_rec.left,
      inp_seqs=inp_seq_rec.right,
      seed=seed,
      read_count=read_count
    }
  }
  output {
    Array[Array[File]] out_seqs = spr.out_seqs
    Array[String] out_SampleID= spr.out_SampleID
  }
}

task subsample_paired_reads {
  String SampleID
  Array[File] inp_seqs
  Float read_count
  String seed=1
  File inp_seq1=inp_seqs[0]
  File inp_seq2=inp_seqs[1]
  ## WDL gets confused in the output section if variable with the same name is
  ## defined as a variable in the outer scope, so we create a copy with unique name
  String inp_spr_SampleID = SampleID
  command {
    seqtk sample -s${seed} ${inp_seq1} ${read_count} > ${inp_spr_SampleID}_1.fastq
    seqtk sample -s${seed} ${inp_seq2} ${read_count} > ${inp_spr_SampleID}_2.fastq
  }
  output {
    Array[File] out_seqs = [inp_spr_SampleID+"_1.fastq",inp_spr_SampleID+"_2.fastq"]
    String out_SampleID = inp_spr_SampleID
  }
}

