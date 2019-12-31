{
    top_dir : std.extVar("top_dir"),
    test_mode : false,
    threads : 4,
    steps : {
        readqc : { 
            config : "fastqc.jsonnet",
            vars : {
                threads : $.threads
            }
        },
        trimming : { 
            config : "bbduk_trim.jsonnet",
            vars : {
                threads : $.threads
            }
        },
        assembly : { 
            config : "bugbuilder.jsonnet",
            vars : {
                threads : $.threads
            }
        },
        variants : { 
            config : "map_variants.jsonnet",
            vars : {
                threads : $.threads
            }
        }
    },
    output_dir : ".",
    input_iterator : {
        input_dir : "FASTQ/*.gz",
        input_seqfile_ext : ".gz",
        forward : "_R1_",
        reverse : "_R2_",
        samp_id_extractor : "^([0-9]+)_",
        sample_id_prefix : "SA"
    },
    makeflow : {
        workflow : "assembly_batch.mkf",
        env : "%s/%s" % [$.top_dir,"env"],
        makeflow_args : "-T torque -B \'-l nodes=1:ppn=%s\'" % $.threads,
        web : false,
        run : false,
        workflow_script : null
    }
}
