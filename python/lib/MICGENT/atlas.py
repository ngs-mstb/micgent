"""Utility functions around generating and transforming input and output files in workflows"""
from MICGENT.py23 import *
from . import util
from . import converters
import os, re, glob
import warnings
from functools import partial

taxa_lineage_prefixes = ['k__', 'p__', 'c__', 'o__', 'f__', 'g__', 's__']

def pick_taxonomy(df):
    df["tx_bin"] = df.checkm_bin_taxonomy_contained.str.split(";")
    df["tx"] = df.taxonomy.str.split(";")
    ##TODO: only use tx_bin when it is deeper than tx; for now, it always appears to be
    tx_bin_def = ~df.tx_bin.isnull()
    df.loc[tx_bin_def,"tx"] = df.loc[tx_bin_def,"tx_bin"]
    #df.loc[df.tx.isnull(),"tx"][:] = tuple(("{}?".format(_) for _ in taxa_lineage_prefixes))

    def _extract_tx_low(x):
        y = [_ for _ in x if not "?" in _] if not isinstance(x, float) else []
        return y[-1] if y else None
    df["tx_low"] = df.tx.apply(_extract_tx_low)
    return df

def extract_func_summary(summ_in,summ_out=None,filter_gene_defined=True):
    import pandas as pd
    summ_df = pd.read_table(summ_in)
    if filter_gene_defined:
        summ_df = summ_df.loc[~summ_df.gene.isnull()]
    summ_df = pick_taxonomy(summ_df)
    if summ_out:
        summ_df.to_csv(summ_out,index=False,sep="\t")
    else:
        return summ_df
    #rep_all = pd.concat(rep_all, axis=0).reset_index(drop=True)

def make_str_extractor(patt):
    def _extract_str(patt,s):
        samp_id = re.findall(patt,s)
        if len(samp_id) < 1:
            return None
        if len(samp_id) > 1 or not samp_id[0].strip():
            raise ValueError("Pattern '{}' for string should extract one "\
                    "or zero non blank matches from '{}', but it extracted this: '{}'".\
                    format(patt,s,samp_id))
        return samp_id[0]
    patt = util.none_from_str(patt)
    if patt:
        if not '(' in patt:
            patt = "({})".format(patt)
        extr = partial(_extract_str,patt)
    else:
        extr = None
    return extr

def extract_func_summary_many(summ_in,summ_out,filter_gene_defined=True,samp_id_extractor=None):
    samp_id_extractor = make_str_extractor(samp_id_extractor)
    for i_file,summ_in_fn in enumerate(util.glob_files(files_globs=summ_in)):
        summ_df = extract_func_summary(summ_in_fn,filter_gene_defined=filter_gene_defined)
        samp_id = os.path.basename(summ_in_fn)
        if samp_id_extractor:
            samp_id = samp_id_extractor(samp_id)
        summ_df["SampleID"] = samp_id
        write_header = False
        mode = "a"
        if i_file == 0:
            write_header = True
            mode = "w"
        summ_df.to_csv(summ_out, index=False, sep="\t", header=write_header, mode=mode)

## import package module and add argh entry points

def _main():
    import argh
    argh.dispatch_commands([
        extract_func_summary,
        extract_func_summary_many
    ])


if __name__ == "__main__":
    _main()