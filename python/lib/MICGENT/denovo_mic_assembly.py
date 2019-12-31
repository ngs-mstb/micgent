"""Pipeline for microbial denovo assembly"""
from __future__ import absolute_import
from __future__ import unicode_literals
from __future__ import print_function
from __future__ import division
from future import standard_library
standard_library.install_aliases()
from builtins import *
from . import util
from . import config
import glob
import os, shutil, re

from MGT.Logging import *
from MGT.BatchMakeflow import makeflow, make_task_sandbox

from . import workflow_util

def extract_bugbuilder_asm_metrics(log_file):
    import re
    def skip_until(inp, rx, n=1):
        i = 0
        for line in inp:
            if re.search(rx, line):
                i += 1
                if i >= n:
                    return line
        return None

    def read_while(inp, rx):
        for line in inp:
            if not re.search(rx, line):
                break
            yield line

    def read_until(inp, rx):
        for line in inp:
            if re.search(rx, line):
                break
            yield line

    tab_spl = r"\s*\|\s*"

    def tab_row(line):
        return re.split(tab_spl, line)[1:-1]

    def tab_rows(inp, n_cols=2):
        for line in read_while(inp, r"^\|"):
            row = tab_row(line)
            if n_cols:
                row = row[:n_cols]
            yield row

    tab_section = r"^\+-"

    def tab_find_and_read(inp, tbl_header, prefix=[], sections_skip=1, n_cols=2):
        skip_until(inp, tbl_header)
        skip_until(inp, tab_section, sections_skip)
        for row in tab_rows(inp, n_cols):
            yield (prefix + row)

    def assembly_statistics(inp, section):
        group = "All Contigs"
        for row in tab_find_and_read(inp, group, prefix=[section, group]):
            yield row
        group = "All Scaffolds"
        for row in tab_find_and_read(inp, group, prefix=[section, group]):
            yield row

    recs = []
    with open(log_file, "r") as inp:
        section = "Read Information"
        group = section
        for row in tab_find_and_read(inp, group, prefix=[section, group], sections_skip=2):
            yield row
        section = "spades assembly statistics"
        skip_until(inp, section)
        for row in assembly_statistics(inp, section):
            yield row
        section = "Final Assembly Statistics"
        skip_until(inp, section)
        for row in assembly_statistics(inp, section):
            yield row
        section = "Annotated features"
        group = "Feature Type"
        for row in tab_find_and_read(inp, group, prefix=[section, group]):
            yield row


def collect_bugbuilder_asm_metrics(log_files, out_file):
    import csv
    header = ("LogFile", "Section", "Group", "Key", "Value")
    with open(out_file, "w") as out:
        wr = csv.writer(out, dialect="excel-tab", lineterminator="\n")
        wr.writerow(header)
        for log_file in util.glob_files(files_globs=log_files):
            for rec in extract_bugbuilder_asm_metrics(log_file):
                wr.writerow([log_file] + rec)


def assemble_de_novo(config_file="assembly_batch.jsonnet"):
    """
    Run de novo assembly  on a set of directories with FASTQ files
    """

    samp_output_dirs = []
    samp_output_files = []

    top_dir = os.getcwd()
    top_dir_readqc_before_trim = os.path.join(top_dir, "readqc_before_trimming")
    top_dir_readqc_after_trim = os.path.join(top_dir, "readqc_after_trimming")
    top_dir_trim = os.path.join(top_dir, "trimming")
    top_dir_asse = os.path.join(top_dir, "assembly")
    top_dir_vari = os.path.join(top_dir, "variants")

    vars_default = dict(top_dir=top_dir)

    batch_conf = config.load_config(config_file, vars_default=vars_default)

    batch_conf["makeflow"].setdefault("wrapper", config.load_config(pkg="MICGENT")["wrapper"])

    with makeflow(**batch_conf["makeflow"]) as mf:

        ci_readqc = config.load_subconfig_interpreter(batch_conf["steps"]["readqc"],
                                                      vars_default=vars_default)
        ci_trim = config.load_subconfig_interpreter(batch_conf["steps"]["trimming"],
                                                    vars_default=vars_default)
        ci_asse = config.load_subconfig_interpreter(batch_conf["steps"]["assembly"],
                                                    vars_default=vars_default)
        ci_vari = config.load_subconfig_interpreter(batch_conf["steps"]["variants"],
                                                    vars_default=vars_default)

        for samp_iid, samp in enumerate(workflow_util.iterate_sequence_run(to_abspath=True,
            **batch_conf["input_iterator"])):

            inp_seqs = samp["seqfiles"]
            samp_id = samp["samp_id"]
            assert not ("_" in samp_id or "." in samp_id), \
                    "Symbols '_.' are not allowed in sample ID: {}".format(samp_id)
            log.info("Creatings tasks for sample {}".format(samp["samp_dir"]))

            id_task = "{}_{}_rq_b_tr".format(samp_id, samp_iid)
            c_task = make_task_sandbox(top_dir_readqc_before_trim,
                                       id_task=id_task)
            out_dir = c_task["cwd"]
            c_task.update(ci_readqc(dict(inp_seqs=inp_seqs,
                                         out_dir=out_dir)))
            mf.task(
                conf=c_task
            )

            id_task = "{}_{}_tr".format(samp_id, samp_iid)
            c_task = make_task_sandbox(top_dir_trim,
                                       id_task=id_task)
            # out_seqs=[ os.path.join(c_task["cwd"],"read_{}.fq.gz".format(_)) for _ in ("1P","1U","2P","2U") ]
            out_seqs = [os.path.join(c_task["cwd"], "read_{}.fq.gz".format(_)) for _ in ("1P", "2P")]
            c_task.update(ci_trim(dict(inp_seqs=inp_seqs,
                                       out_seqs=out_seqs)))
            mf.task(
                conf=c_task
            )

            id_task = "{}_{}_rq_a_tr".format(samp_id, samp_iid)
            c_task = make_task_sandbox(top_dir_readqc_after_trim,
                                       id_task=id_task)
            out_dir = c_task["cwd"]
            c_task.update(ci_readqc(dict(inp_seqs=out_seqs,
                                         out_dir=out_dir)))
            mf.task(
                conf=c_task
            )

            id_task = "{}_{}_as".format(samp_id, samp_iid)
            c_task = make_task_sandbox(top_dir_asse,
                                       id_task=id_task)
            inp_seqs = [out_seqs[0], out_seqs[1]]
            out_dir = os.path.join(c_task["cwd"], "out")
            c_task.update(ci_asse(dict(inp_seqs=inp_seqs,
                                       out_dir=out_dir,
                                       locustag=samp_id)))
            mf.task(
                conf=c_task
            )

            id_task = "{}_{}_va".format(samp_id, samp_iid)
            c_task = make_task_sandbox(top_dir_vari,
                                       id_task=id_task)
            inp_seqs = [out_seqs[0], out_seqs[1]]
            out_bam = os.path.join(c_task["cwd"], "final.bam")
            out_bcf = os.path.join(c_task["cwd"], "final.bcf")
            c_task.update(ci_vari(dict(inp_seqs=inp_seqs,
                                       out_bam=out_bam,
                                       out_bcf=out_bcf)))
            mf.task(
                conf=c_task
            )

            if batch_conf["test_mode"]:
                break


def post_assembly_fix(asm_dir, samp_id=""):
    """Apply a few post hoc fixes to one assembly directory.
    This currently turns scaffold names into prefix_samplid_scaffold_number names.
    and generates Genabk files from EMBL. Fixes from here should be eventualy
    migrated into the assembly pipeline.
    @param asm_dir directory with output assembly files (such as BugBuilder `out` dir)
    @param samp_id Should look like "STRAINID". STRAINID should be globally unique
    but short.
    """
    from subprocess import check_call
    if not samp_id:
        samp_id = "S"
    out_dir = asm_dir
    target_files = []
    for ext in (".fasta", ".embl", ".agp"):
        target_files += list(glob.glob(os.path.join(out_dir, "*" + ext)))
    check_call(["perl", "-p", "-i", "-e", 's/scaffold_/{}_/g'.format(samp_id)] + \
               target_files)
    check_call(["perl", "-p", "-i", "-e", 's/contig_/{}_ctg_/g'.format(samp_id)] + \
               target_files)
    embl_file = os.path.join(out_dir, "scaffolds.embl")
    gb_file = os.path.join(out_dir, "scaffolds.gb")
    check_call(["seqret", "-sequence", "embl::{}".format(embl_file), "-outseq", gb_file,
                "-osformat2", "genbank", "-feature"])
    gff_file = os.path.join(out_dir, "scaffolds.gff")
    check_call(["seqret", "-sequence", "embl::{}".format(embl_file), "-outseq", gff_file,
                "-osformat2", "gff", "-feature"])


def post_assembly_fix_many(asm_dirs, samp_id_prefix=""):
    """Apply post_assembly_fix to the collection of assembly directories.
    @param glob for directories, one per BugBuilder run for a single genome (each has subdir `out`)
    @param prefix Should look like "EXPERIMENT_". EXPERIMENT should be globally unique
    but short.
    """
    import glob
    for asm_dir in util.glob_files(files_globs=asm_dirs):
        samp_id = os.path.basename(asm_dir).split("_")[0]
        out_dir = os.path.join(asm_dir, "out")
        post_assembly_fix(asm_dir=out_dir, samp_id="{}{}".format(samp_id_prefix,samp_id))

def export_assembly(asm_dirs, out_dir, export_methods="hsc"):
    """Export assembly output files from many directories into one.
    @param glob for directories with the output files, one per genome
    @param out_dir output directory
    @param export_methods Try in the given order hard linking, soft 
    linking or copying.
    @note The first part of the FASTA ID (before '_') will be used to make
    the file name (before the .extension)
    """
    import glob
    util.makedir(out_dir)
    for asm_dir in util.glob_files(files_globs=asm_dirs):
        scaf_file = os.path.join(asm_dir,"scaffolds.fasta")
        with open(scaf_file,"r") as inp:
            defline = inp.readline()
            samp_id = defline.split(">")[1].split("_")[0]
        for fn in glob.glob(os.path.join(asm_dir,'*')):
            if os.path.isfile(fn):
                fn_base = os.path.basename(fn)
                fn_base_root,fn_base_ext = fn_base.rsplit(".",1)
                if fn_base_root == "scaffolds":
                    fn_base_root = samp_id
                else:
                    fn_base_root = samp_id+"."+fn_base_root
                fn_out_base = fn_base_root+"."+fn_base_ext
                fn_out = os.path.join(out_dir,fn_out_base)
                util.link_or_copy(os.path.abspath(fn),fn_out)
                

def re_export_assembly(asm_files, 
        out_dir, 
        id_file, 
        id_column="organism",
        negate=False,
        id_column_prefix_add="",
        export_methods="hsc"):
    """Re-export a subset of assembly output files from one directory to another.
    @param asm_files glob for asm_files, can be spread across multiple directories,
    each with output files for many assemblies (names must be globally unique).
    Such directory can be created with export_assembly method.
    @param id_file CSV file with IDs to select
    @param column name in CSV file to use for IDs selection
    @negate If True, export everything but what is in id_file [False]
    @param id_column_prefix_add Add this prefix to value in ID column
    @param out_dir output directory
    @param export_methods Try in the given order hard linking, soft 
    linking or copying.
    """
    import glob
    import pandas as pd
    ids = set([ "{}{}".format(id_column_prefix_add,_) \
            for _ in pd.read_csv(id_file)[id_column]])
    util.makedir(out_dir)
    for asm_file in util.glob_files(files_globs=asm_files):
        fn_base = os.path.basename(asm_file)
        seq_id = fn_base.split(".")[0]
        if (not negate and seq_id in ids) or (negate and seq_id not in ids):
            fn_out = os.path.join(out_dir,fn_base)
            util.link_or_copy(os.path.abspath(asm_file),fn_out)
                


## import package module and add argh entry points

def _main():
    import argh
    argh.dispatch_commands([
        assemble_de_novo,
        collect_bugbuilder_asm_metrics,
        export_assembly,
        re_export_assembly,
        post_assembly_fix,
        post_assembly_fix_many
    ])


if __name__ == "__main__":
    _main()
