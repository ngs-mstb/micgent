"""Wrappers to run sequence QC, trimming and assembly is distributed mode.
Uses Makeflow.
"""
from __future__ import absolute_import
from __future__ import unicode_literals
from __future__ import print_function
from __future__ import division

from future import standard_library
standard_library.install_aliases()
from builtins import *
from . import config
import glob
import os, shutil, re

from MGT.Logging import *
from MGT.BatchMakeflow import makeflow, make_task_sandbox

def trim_galore(input_dir,
        stringency=1,
        length=75,
        no_fastqc=False,
        trim_galore_args="",
        fastqc_args="",
        output_dir=".",
        wrapper=None,
        input_seqfile_ext=".fastq.gz",
        forward="_R1_001",
        reverse="_R2_001",
        makeflow_args="",
        web=False,
        run_no=False):
    """
    Run trim_galore on a set of directories with FASTQ files
    """

    if not wrapper:
        wrapper = config.load_config(pkg="MICGENT")["wrapper"]

    ## --other argument string to srst2.py need involved quoting of spaces
    ## in order to pass through Makeflow shell first, and then through
    ## ArgParse splitting untouched
    if fastqc_args.strip():
        fastqc_args = [fastqc_args.strip()]
    else:
        fastqc_args = []
    fastqc_args = r"\\ ".join(fastqc_args)
    if fastqc_args:
        fastqc_args = "--fastqc_args '{}'".format(fastqc_args)

    if no_fastqc:
        fastqc = ""
        fastqc_args = ""
    else:
        fastqc = "--fastqc"

    trim_galore_args += " --stringency {} --length {}".format(stringency,length)

    samp_iid = 0
    samp_output_dirs = []
    samp_output_files = []
    
    workflow = "trim_galore.mkf"

    run_pref = "trim_galore"

    with makeflow(
        workflow=workflow,
        wrapper=wrapper,
        makeflow_args=makeflow_args,
        workflow_script=None,
        web=web,
        run=not run_no
        ) as mf:

        ##TODO: add first step that indexes the reference, as SRST2 manual advises 
        for samp_dir in glob.glob(input_dir):
            if os.path.isdir(samp_dir):

                samp_trim_galore_args = trim_galore_args

                samp_seqfiles = glob.glob(os.path.join(samp_dir,"*"+input_seqfile_ext))
                assert len(samp_seqfiles), \
                        "No sequence files found for subdirectory {}".format(samp_dir)
                if len(samp_seqfiles) > 1:
                    ## reorder files in forward, reverse pattern order
                    assert len(samp_seqfiles) == 2, "Expecting one pair of FASTQ files per sample directory"
                    for ff in (samp_seqfiles,tuple(reversed(samp_seqfiles))):
                        if re.search(forward,ff[0]):
                            assert re.search(reverse,ff[1]),"Forward read file name pattern matched, but reverse pattern did not: {}".format(ff)
                        samp_seqfiles = ff
                        break
                    else:
                        raise ValueError("Did not find forward and reverse file name patterns: {}".format(samp_seqfiles))
                    samp_trim_galore_args += " --paired"

                samp_output_pref = "{}_{}".format(run_pref,samp_iid)
                samp_output_dir = os.path.join(output_dir,samp_output_pref)
                shutil.rmtree(samp_output_dir,ignore_errors=True)
                os.makedirs(samp_output_dir)
                status_file = samp_output_pref+".workflow_status"
            
                cmd = """trim_galore \
                --output_dir {output_dir} \
                {fastqc} \
                {fastqc_args} \
                {trim_galore_args} \
                {seqfiles}""".format(seqfiles=" ".join(samp_seqfiles),
                        output_dir=samp_output_dir,
                        trim_galore_args=samp_trim_galore_args,
                        fastqc=fastqc,
                        fastqc_args=fastqc_args)
                
                mf.task(
                        cmd=cmd,
                        inputs=samp_seqfiles,
                        status_file=status_file
                        )

                samp_output_dirs.append(samp_output_dir)
                samp_output_files.append(status_file)
                samp_iid += 1
    
        
        if False:
            cmd = """echo \
                  {samp_output_dirs}""".format(samp_output_dirs=" ".join(samp_output_dirs))
                
            mf.task(
                    cmd=cmd,
                    targets=[],
                    inputs=samp_output_files
                    )

## import package module and add argh entry points

def _main():
    import argh
    argh.dispatch_commands([
        trim_galore
    ])


if __name__ == "__main__":
    _main()
