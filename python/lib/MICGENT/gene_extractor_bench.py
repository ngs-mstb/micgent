#!/usr/bin/env python
"""Some benchmarking and testing helper methods for the gene extractor"""

from __future__ import print_function
from __future__ import unicode_literals
from __future__ import division
from __future__ import absolute_import
from future import standard_library
standard_library.install_aliases()
from builtins import *
from . import util
from . import workflow_util
from . import sig

import pandas as pd
import numpy as np

import os
import itertools

from MGT.Logging import *


class DeterminismError(Exception):
    """Base class for exceptions when checks of determinism fail."""
    def __init__(self,context,metrics):
        self.context = context
        self.metrics = metrics
    def __str__(self):
        return "Determinism checks failed for {}, with metrics {}".format(self.context,self.metrics)

def replicate_sample_man(man_inp,n_inp=1,n_repl=10):
    """Replicate first n_inp samples from the manifest into multiple copies with different SampleID fields"""
    man = pd.read_table(man_inp)
    if n_inp>0:
        man = man.head(n_inp)
    man_out = None
    for i in range(1,n_repl+1):
        man_repl = man.copy()
        man_repl.SampleID = man_repl.SampleID.apply(lambda x: "{}_{}".format(x,i))
        if man_out is None:
            man_out = man_repl
        else:
            man_out = man_out.append(man_repl,ignore_index=True)
    return man_out

def expand_samp_arch_dirs(top_dir,arch_glob="*.ariba.debug.tar",skip_safety_check=False):
    for gl in arch_glob.split(","):
        for arch_file in util.glob_files(os.path.join(top_dir,gl)):
            arch_dir = arch_file+".dir"
            util.tar_to_dir(arch_file,arch_dir,skip_safety_check=skip_safety_check,overwrite=True)

def compute_sigs_run(run_dir,arch_glob="*.ariba.debug.tar",out_csv="sig_run.txt"):
    ## skipping safety check to allow symlinks - the archives in this module should
    ## only come from the internal runs
    expand_samp_arch_dirs(run_dir,arch_glob=arch_glob,skip_safety_check=True)
    sig.dir_sig(run_dir,with_dir_name=True,no_recurse=False,out_csv=out_csv)

def run_out_files_to_check(sor_pack_dir):
    """Return a list of output file paths to check, relative to unspecified top dir
    
    :param sor_pack_dir: path to the directory of the expanded sor pack, relative to
    the top output directory.
    
    No check is made that files exist as this point.
    """
    manifest_files = ["manifest_out.tsv",
    "web/manifest_out_all.tsv",
    os.path.join(sor_pack_dir,"manifest.tsv")]
    fasta_files = ["seq_out.fasta",
    "web/seq_out_all.fasta",
    os.path.join(sor_pack_dir,"seq.fasta")]
    return (manifest_files,fasta_files)

def extract_non_repl_id(x):
    if pd.notnull(x):
        parts = x.rsplit("_",2)
        x = "_".join((parts[0],parts[-1]))
    return x

def check_table_replicate_equality(man,no_assert=False,assert_msg=""):
    man = man.copy()
    ## SampID structure is assumed to be created by replicate_sample_man
    man["SampleID_Base"] = man.SampleID.str.rsplit("_",1).str[0]
    ## pandas.groupby silently drops NA key values, so we replace NA SeqID with
    ## a sample-specific string value in order to get them counted and checked for
    ## agreement between replicates
    null_SeqID = man.SeqID.isnull()
    if null_SeqID.any():
        #import pdb; pdb.set_trace()
        man.loc[null_SeqID,"SeqID"] = man.SampleID[null_SeqID].apply(lambda x: "{}_NA".format(x))
    ## Derive SeqID that would be expected from the original non-replicated
    ## SampleID
    man["SeqID_Base"] = man.SeqID.apply(extract_non_repl_id) 
    ## In a copy of man, replace SampleID and SeqID with base values,
    ## and then test that there is a single value for every remaining field
    ## for all records (replicates) with a given value of SeqID (which is now base SeqID)
    man_rep = man.copy()
    man_rep.SampleID = man_rep.SampleID_Base
    man_rep.SeqID = man_rep.SeqID_Base
    if "SeqID_SOR" in man_rep:
        man_rep.SeqID_SOR = man_rep.SeqID_SOR.apply(extract_non_repl_id)
    ## all values for  a given base SeqID are either equal and non-null or just null
    is_eq_df = man_rep.groupby("SeqID").agg(lambda x:(np.unique(x).size==1 and pd.notnull(x).all()) or pd.isnull(x).all())
    is_eq_flds = is_eq_df.all()
    non_eq_fld_recs = is_eq_df[is_eq_df.sum(axis=1)<is_eq_df.shape[1]]
    ## this is the final bool for field equality check
    is_eq_all_flds = is_eq_flds.all()
    ## Now check that all replicates of a given sample generated equal number of sequence records
    uni_seq_counts = man.groupby(["SampleID_Base","SampleID"]).\
        SeqID.count().groupby("SampleID_Base").unique()
    ## this is the final bool for sequence count per sample per replicate equality check
    all_eq_seq_count = (uni_seq_counts.apply(len) == 1).all()
    n_samp = len(man.SampleID.unique())
    n_samp_base = len(man.SampleID_Base.unique())
    res = dict(is_eq_all_flds=is_eq_all_flds,
        is_eq_flds=is_eq_flds,
        is_eq_df=is_eq_df,
        non_eq_fld_recs=non_eq_fld_recs,
        uni_seq_counts=uni_seq_counts,
        all_eq_seq_count=all_eq_seq_count,
        n_samp=n_samp,
        n_samp_base=n_samp_base,
        n_rec=man.shape[0])
    if not no_assert:
        if not res["is_eq_all_flds"] or not res["all_eq_seq_count"]:
            raise DeterminismError(assert_msg,res)
    if man.SeqID.isnull().any():
        import pdb; pdb.set_trace()
    return res

def check_manifest_replicate_equality(man_file,no_assert=False):
    man = workflow_util.load_manifest(man_file)
    return check_table_replicate_equality(man,no_assert=no_assert,assert_msg=man_file)

def check_fasta_replicate_equality(fasta_file,out_sig_csv=None,strip_sfx=None,no_assert=False):
    if out_sig_csv is None:
        out_sig_csv = fasta_file+".sig.txt"
    sig.bioseq_sig(fasta_file,out_csv=out_sig_csv)
    sg = pd.read_table(out_sig_csv,header=None,names=["SeqID","Sig"],dtype={"SeqID":str})
    if strip_sfx:
        assert sg.SeqID.str.endswith(strip_sfx).all()
        sg.SeqID = sg.SeqID.str.rstrip(strip_sfx)
    sg["SampleID"] = sg.SeqID.str.rsplit("_",1).str[0]
    return check_table_replicate_equality(sg,no_assert=no_assert,assert_msg=fasta_file) 

def unzip_sor_pack(sor_pack_file="sor_pack.tgz",sor_pack_dir="sor_pack"):
    util.tar_to_dir(sor_pack_file,sor_pack_dir,skip_safety_check=True,overwrite=True)
    return sor_pack_dir


def check_within_run_replicate_equality(run_out_dir,no_assert=False):
    with util.chdir(run_out_dir,create=False,to_abs=True):   
        sor_pack_dir = unzip_sor_pack()
        (manifest_files,fasta_files) = run_out_files_to_check(sor_pack_dir)
        out_man = {}
        for manifest_file in manifest_files:
            out_man[manifest_file] = check_manifest_replicate_equality(manifest_file,no_assert=no_assert)
        out_fasta = {}
        for fasta_file in fasta_files:
            out_fasta[fasta_file] = check_fasta_replicate_equality(fasta_file,no_assert=no_assert)
        return dict(manifest=out_man,
            fasta=out_fasta)

def compute_run_sigs(run_out_dir,out_csv,use_existing=False):
    if use_existing and os.path.exists(out_csv):
        df = pd.read_table(out_csv)
    else:
        recs = []
        with util.chdir(run_out_dir,create=False,to_abs=True):   
            sor_pack_dir = unzip_sor_pack()
            (manifest_files,fasta_files) = run_out_files_to_check(sor_pack_dir)
            for fn in itertools.chain(manifest_files,fasta_files):
                rec = None
                if not os.path.isfile(fn):
                    rec = (fn,None)
                else:
                    rec = (fn,sig.file_sig(fn))
                recs.append(rec)
        df = pd.DataFrame.from_records(recs,columns=["File","Signature"])
        df.to_csv(out_csv,index=False)
    return df

def check_runs_for_equality(run_out_dirs,no_assert=False,use_existing_sigs=False):
    runs_sigs = [ ]
    for run_out_dir in util.glob_files(run_out_dirs,no_glob_match="raise"):
        runs_sigs.append((run_out_dir,compute_run_sigs(run_out_dir,
            os.path.join(run_out_dir,"out_sigs.txt"),
            use_existing=use_existing_sigs)))

    assert len(runs_sigs) > 1, "Need at least two run output directories to compare"
    run_one = runs_sigs[0]
    runs_rest = runs_sigs[1:]
    cmps = []
    for run_two in runs_rest:
        pair_eq = run_one[1].equals(run_two[1])
        if not pair_eq and not no_assert:
            raise DeterminismError((run_one[0],run_two[0]),
                (run_one[1],run_two[1]))
        cmps.append((run_one[0],run_two[0],pair_eq))
    return cmps

## import package module and add argh entry points

def _main():
    import argh
    argh.dispatch_commands([
        replicate_sample_man,
        compute_sigs_run,
        check_manifest_replicate_equality,
        check_fasta_replicate_equality,
        check_within_run_replicate_equality,
        compute_run_sigs,
        check_runs_for_equality        
    ])


if __name__ == "__main__":
    _main()
