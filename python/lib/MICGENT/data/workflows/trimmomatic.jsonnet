{
    top_dir  : std.extVar("top_dir"),
    adapter_files : [ "%s/%s" % [self.top_dir,f] for f in ["TruSeq3-PE-2.fa","TruSeq.AdInd.fa"] ],
    threads : std.extVar("threads"),
    inp_seqs : std.split(std.extVar("inp_seqs"),","),
    out_seqs : std.split(std.extVar("out_seqs"),","),

    cmd : |||
        trimmomatic PE -threads %s \
        -trimlog trim.log %s %s \
        ILLUMINACLIP:%s:2:30:10:1:true \
        ILLUMINACLIP:%s:0:11:6:1:true
    ||| % [self.threads,
    std.join(" ",self.inp_seqs),
    std.join(" ",self.out_seqs),
    self.adapter_files[0],
    self.adapter_files[1]],
    // LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:20

    inputs : self.inp_seqs + self.adapter_files,
    targets : self.out_seqs
}

