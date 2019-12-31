{
    top_dir  : std.extVar("top_dir"),
    inp_seqs : std.split(std.extVar("inp_seqs"),","),
    inp_seq1 : self.inp_seqs[0],
    inp_seq2 : self.inp_seqs[1],
    out_dir  : std.extVar("out_dir"),
    locustag : std.extVar("locustag"),
    threads : std.extVar("threads"),
    genus : "Staphylococcus",
    species : "aureus",
    reference : "%s/%s" % [self.top_dir,"NC_007795.fna"],

    tmp_dir : self.out_dir + ".tmp",

    cmd : |||
        BugBuilder --no-gap-fill \
        --mode draft \
        --assembler-args=\'-t %(threads)s --careful -k 21,33,55,77,99,127\' \
        --out-dir %(out_dir)s \
        --tmp-dir %(tmp_dir)s \
        --scaffolder ragout \
        --genus %(genus)s --species %(species)s \
        --locustag %(locustag)s \
        --threads %(threads)s \
        --reference %(reference)s \
        --fastq1 %(inp_seq1)s \
        --fastq2 %(inp_seq2)s \
        --platform illumina \
        --assembler spades \
        --min-contig 300 \
        --no-fastqc \
        --no-reorder-scaffolds \
        --no-split-origin && \
        python -m MICGENT.denovo_mic_assembly post-assembly-fix --samp-id %(locustag)s %(out_dir)s
    ||| % self,
    cwd_remove_before : true,
    inputs : self.inp_seqs + [ self.reference ],
    targets : [ self.out_dir + "/" + f for f in ["scaffolds.agp","scaffolds.fasta",
                "circleator.png","assembly.bnk","scaffolds.embl",
                "scaffolds.gb","scaffolds.gff"] ]
}

