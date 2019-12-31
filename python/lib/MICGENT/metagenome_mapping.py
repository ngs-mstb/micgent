"""Pipeline for mapping metagenomic WGS reads to reference"""
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


def metagenome_map(config_file="metagenome_map_batch.jsonnet"):
    """
    Map metagenomes from a set of directories with FASTQ files
    """

    samp_output_dirs = []
    samp_output_files = []

    top_dir = os.getcwd()
    top_dir_readqc_before_trim = os.path.join(top_dir, "readqc_before_trimming")
    top_dir_readqc_after_trim = os.path.join(top_dir, "readqc_after_trimming")
    top_dir_trim = os.path.join(top_dir, "trimming")
    top_dir_map = os.path.join(top_dir, "map")

    vars_default = dict(top_dir=top_dir)

    batch_conf = config.load_config(config_file, vars_default=vars_default)

    batch_conf["makeflow"].setdefault("wrapper", config.load_config(pkg="MICGENT")["wrapper"])

    with makeflow(**batch_conf["makeflow"]) as mf:

        ci_readqc = config.load_subconfig_interpreter(batch_conf["steps"]["readqc"],
                                                      vars_default=vars_default)
        ci_trim = config.load_subconfig_interpreter(batch_conf["steps"]["trimming"],
                                                    vars_default=vars_default)
        ci_map = config.load_subconfig_interpreter(batch_conf["steps"]["map"],
                                                   vars_default=vars_default)

        for samp_iid, samp in enumerate(iterate_sequence_run(**batch_conf["input_iterator"])):

            inp_seqs = samp["seqfiles"]
            samp_base = os.path.basename(samp["samp_dir"])
            log.info("Creatings tasks for sample {}".format(samp["samp_dir"]))

            id_task = "{}_{}_rq_b_tr".format(samp_base, samp_iid)
            c_task = make_task_sandbox(top_dir_readqc_before_trim,
                                       id_task=id_task)
            out_dir = c_task["cwd"]
            c_task.update(ci_readqc(dict(inp_seqs=inp_seqs,
                                         out_dir=out_dir)))
            mf.task(
                conf=c_task
            )

            id_task = "{}_{}_tr".format(samp_base, samp_iid)
            c_task = make_task_sandbox(top_dir_trim,
                                       id_task=id_task)
            # out_seqs=[ os.path.join(c_task["cwd"],"read_{}.fq.gz".format(_)) for _ in ("1P","1U","2P","2U") ]
            out_seqs = [os.path.join(c_task["cwd"], "read_{}.fq.gz".format(_)) for _ in ("1P", "2P")]
            c_task.update(ci_trim(dict(inp_seqs=inp_seqs,
                                       out_seqs=out_seqs)))
            mf.task(
                conf=c_task
            )

            id_task = "{}_{}_rq_a_tr".format(samp_base, samp_iid)
            c_task = make_task_sandbox(top_dir_readqc_after_trim,
                                       id_task=id_task)
            out_dir = c_task["cwd"]
            c_task.update(ci_readqc(dict(inp_seqs=out_seqs,
                                         out_dir=out_dir)))
            mf.task(
                conf=c_task
            )

            id_task = "{}_{}_ma".format(samp_base, samp_iid)
            c_task = make_task_sandbox(top_dir_map,
                                       id_task=id_task)
            inp_seqs = [out_seqs[0], out_seqs[1]]
            out_bam = os.path.join(c_task["cwd"], "final.bam")
            out_bcf = os.path.join(c_task["cwd"], "final.bcf")
            c_task.update(ci_map(dict(inp_seqs=inp_seqs,
                                      out_bam=out_bam,
                                      out_bcf=out_bcf)))
            mf.task(
                conf=c_task
            )

            if batch_conf["test_mode"]:
                break


## import package module and add argh entry points

def _main():
    import argh
    argh.dispatch_commands([
        metagenome_map
    ])


if __name__ == "__main__":
    _main()
