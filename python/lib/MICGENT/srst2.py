"""Wrappers to run SRST2 is distributed mode.
SRST2 (https://github.com/katholt/srst2) annotates MLST, resistance genes or anything that is based on nucleotide
reference database from raw reads.
These wrappers automate running SRST2 in parallel on multiple input datasets and
combining the results using Makeflow.
"""
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
import os

from MGT.Logging import *
from MGT.BatchMakeflow import makeflow

def srst_output_pref_to_file_name(pref,file_type="mlst_results",resolve_glob=False):
    ret = None
    if file_type == "mlst_results":
        ret = "{}__mlst__*__results.txt".format(pref)
    elif file_type == "compiled_results":
        ret = "{}__compiledResults.txt".format(pref)
    else:
        raise ValueError("Unknown file_type: {}".format(file_type))
    if resolve_glob:
        ret = glob.glob(ret)
        assert len(ret) == 1, "One and only one result file is expected, but instead we found: {}".format(ret)
        ret = ret[0]
    return ret

def mlst(input_dir,
        srst2_args="",
        srst2_other_args="",
        wrapper=None,
        input_seqfile_ext=".fastq.gz",
        forward="_R1_001",
        reverse="_R2_001",
        no_discordant=False,
        no_contained=False,
        mlst_max_mismatch=None,
        makeflow_args="",
        web=False,
        run_no=False):
    """
    python -m MICGENT.srst2.py mlst --wrapper ./wrapper \
            --srst2-args "--mlst_db Klebsiella_pneumoniae.fasta --mlst_definitions kpneumoniae.txt --mlst_delimiter '_'" \
            --srst2-other-args "\\\\--no-discordant\\\\ --no-contain" \
            --makeflow-args "-T torque" \
            "$(pwd)/samples/*"
    Note the quad quoting of spaces with backslash between elements of the --srst2-other-args argument. By default, we add
    --no-discordant and --no-contain to the SRST2 --other argument to mitigate cases when partial matches within reads clearly
    belonging to other genes get mapped and downgrade the correct allele, resulting in the incorrect assignment.
    and ./wrapper looks like (this is very specific to the cluster environment - the wrapper should make sure that srst executable
    and all its dependencies are available to the cluster job when it runs):
    ```

    ```
    """

    if not wrapper:
        wrapper = config.load_config(pkg="MICGENT")["wrapper"]

    ## --other argument string to srst2.py need involved quoting of spaces
    ## in order to pass through Makeflow shell first, and then through
    ## ArgParse splitting untouched
    if srst2_other_args.strip():
        srst2_other_args = [srst2_other_args.strip()]
    else:
        srst2_other_args = []
    if no_discordant:
        srst2_other_args.append("--no-discordant")
    if no_contained:
        srst2_other_args.append("--no-contain")
    srst2_other_args = r"\\ ".join(srst2_other_args)
    if srst2_other_args:
        srst2_other_args = "--other '{}'".format(srst2_other_args)

    if mlst_max_mismatch is not None:
        mlst_max_mismatch = "--mlst_max_mismatch {}".format(int(mlst_max_mismatch))
    else:
        mlst_max_mismatch = ""
    samp_iid = 0
    samp_output_prefs = []
    srst_output_files = []
    
    workflow = "mlst.mkf"

    run_pref = "micgent_samp"

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

                samp_seqfiles = glob.glob(os.path.join(samp_dir,"*"+input_seqfile_ext))
                assert len(samp_seqfiles), \
                        "No sequence files found for subdirectory {}".format(samp_dir)

                samp_output_pref = "{}_{}".format(run_pref,samp_iid)
                samp_output_file = srst_output_pref_to_file_name(samp_output_pref,file_type="mlst_results")
                status_file = samp_output_pref+".workflow_status"
            
                ## SRST2 does not exit with error code if it detects an error condition,
                ## so we need to check for file existence.
                ## Because this will be passed to wrapper script,
                ## we need to make it into a single command as bash
                ## argument (brackets will not work).
                cmd = """srst2 \
                --input_pe {seqfiles} \
                --output {output_pref} \
                --forward {forward} \
                --reverse {reverse} \
                {mlst_max_mismatch} \
                {srst2_args} {srst2_other_args} && test -f {output_file}""".format(seqfiles=" ".join(samp_seqfiles),
                        output_pref=samp_output_pref,
                        output_file=samp_output_file,
                        srst2_args=srst2_args,
                        srst2_other_args=srst2_other_args,
                        forward=forward,
                        reverse=reverse,
                        mlst_max_mismatch=mlst_max_mismatch)
                
                mf.task(
                        cmd=cmd,
                        inputs=samp_seqfiles,
                        status_file=status_file
                        )

                samp_output_prefs.append(samp_output_pref)
                srst_output_files.append(status_file)
                samp_iid += 1
    
        samp_output_globs = " ".join([
            srst_output_pref_to_file_name(x,file_type="mlst_results") for x in samp_output_prefs
            ])
        cmd = """srst2 \
              --prev_output {samp_output_globs} \
              --output {run_pref}""".format(samp_output_globs=samp_output_globs,
                      run_pref=run_pref)
            
        mf.task(
                cmd=cmd,
                targets=[srst_output_pref_to_file_name(run_pref,file_type="compiled_results")],
                inputs=srst_output_files
                )

## import package module and add argh entry points

def _main():
    import argh
    argh.dispatch_commands([
        mlst
    ])


if __name__ == "__main__":
    _main()
