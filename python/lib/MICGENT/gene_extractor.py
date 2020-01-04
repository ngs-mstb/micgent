#!/usr/bin/env python

from __future__ import print_function
from __future__ import unicode_literals
from __future__ import division
from __future__ import absolute_import
from future import standard_library
standard_library.install_aliases()
from builtins import *
from . import util
from . import seq_util
from . import num_util
from . import resources
from . import workflow_util
from . import yaml_util
from . import cwl_runner
from . import arg_parsing

import argh
import os
import sys
import shutil
import subprocess
import re
import time
from collections import OrderedDict

import pandas as pd
import numpy as np
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet.IUPAC import ambiguous_dna
from Bio.Seq import Seq

from MGT.Logging import *


def parse_seqmap(seqmap):
    ## Note that match and gene coords have very different meaning in revcomp case - see
    ## extract_padded_* methods below
    seqmap = pd.concat([seqmap,seqmap.seq_id.str.\
        extract((r'(?P<contig>cluster(?:_[0-9]+)*\.l[0-9]+\.c[0-9]+\.ctg\.[0-9]+)'
            r'(?:\.(?P<start>[0-9]+)[-.]'
            r'(?P<end>[0-9]+)\.'
            r'(?P<strand>[+]*))*'),expand=True)
        ],axis=1)
    seqmap.loc[seqmap.start.isnull(),"start"] = "0"
    seqmap["start"] = pd.to_numeric(seqmap.start)
    seqmap.loc[seqmap.end.isnull(),"end"] = "0"
    seqmap["end"] = pd.to_numeric(seqmap.end)
    seqmap["priority"] = 3
    seqmap.loc[seqmap.seq_type=="gene","priority"] = 1
    seqmap.loc[seqmap.seq_type=="match","priority"] = 2
    return seqmap

def load_assemblies_as_df(assembled_genes,assembled_seqs,assemblies,alphabet=ambiguous_dna):
    dfs = []
    df = seq_util.fasta_to_df(assembled_genes,alphabet=alphabet)
    df["seq_type"] = "gene"
    dfs.append(df)
    df = seq_util.fasta_to_df(assembled_seqs,alphabet=alphabet)
    df["seq_type"] = "match"
    dfs.append(df)
    df = seq_util.fasta_to_df(assemblies,alphabet=alphabet)
    df["seq_type"] = "asm"
    dfs.append(df)
    return pd.concat(dfs,axis=0)

def load_seqmap_with_seqs(seq_map,assembled_genes,assembled_seqs,assemblies):
    seqs = load_assemblies_as_df(assembled_genes=assembled_genes,
        assembled_seqs=assembled_seqs,
        assemblies=assemblies)
    assert seqs.shape[0]==0 or seqs.id.is_unique, "Sequence ID must be unique across all sequence files"
    seqmap = parse_seqmap(pd.read_table(seq_map))
    assert seqmap.shape[0]==0 or seqmap.seq_id.is_unique, "Sequence ID must be unique in seqmap table" 
    if seqmap.shape[0] > 0:
        mrg = pd.merge(seqmap,seqs,
            left_on=["seq_id"],
            right_on=["id"],
            how="inner",
            suffixes=["","_seq"])
        mrg.reset_index(drop=True,inplace=True)
    else:
        mrg = seqmap
    assert num_util.series_equal(mrg.seq_id,seqmap.seq_id), "1-to-1 relation is expected for seqmap and sequences"
    return mrg

def prune_seqmap(seqmap):
    seqmap["gene_or_match"] = (seqmap.seq_type != "asm")
    seqmap["seq_len_orig"] = seqmap.seq.str.len()
    seqmap["pos_strand"] = (seqmap.strand == "+")
    ## The next two statements will leave one match or gene per contig, with
    ## a preference for genes, and then for longer sequences and positive strand
    seqmap = seqmap.sort_values(["priority","seq_len_orig","pos_strand"],
        ascending=[True,False,False])
    return seqmap.drop_duplicates(["ref_name","contig","gene_or_match"],keep="first")

def extract_padded_match(rec,pad):
    ## python clips out of bound high index, but negative low reverts 
    ## the coord
    ## from ariba.assembly_compare - it first extracts, then revcomps,
    ## and returns original coords

    start = rec.start - 1
    end = rec.end
    seq = rec.seq_asm[max(0,start-pad):(end+pad)]
    if rec.strand == "-":
        seq = str(Seq(seq, ambiguous_dna).reverse_complement())
    return seq

def extract_padded_gene(rec,pad):
    ## python clips out of bound high index, but negative low reverts 
    ## the coord
    ## from ariba.assembly_compare - it first revcomps, then extracts,
    ## and returns the coords on the revcomped sequence but with start and end swapped

    if rec.start <= rec.end:
        start = rec.start - 1
        end = rec.end
        seq_asm = rec.seq_asm
    else:
        start = rec.end - 1
        end = rec.start
        seq_asm = str(Seq(rec.seq_asm, ambiguous_dna).reverse_complement())

    return seq_asm[max(0,start-pad):(end+pad)]

def seqmap_extend_matches(seqmap,pad_assembled,pad_gene):
    mrg = pd.merge(seqmap,
        seqmap.loc[seqmap.seq_type=="asm",["contig","seq"]],
        on=["contig"],
        how="inner",
        suffixes=["","_asm"])
    mrg.reset_index(drop=True,inplace=True)
    assert num_util.series_equal(mrg.seq_id,seqmap.seq_id), "Expected one asm seq_type for every other seq_type"
    seqmap = mrg
    seqmap["seq_pad"] = seqmap.seq
    if pad_assembled > 0:
        seqmap.loc[seqmap.seq_type=="match","seq_pad"] = seqmap[seqmap.seq_type=="match"].\
        apply(lambda x: extract_padded_match(x,pad=pad_assembled),axis=1)
    if pad_gene > 0:
        seqmap.loc[seqmap.seq_type=="gene","seq_pad"] = seqmap[seqmap.seq_type=="gene"].\
        apply(lambda x: extract_padded_gene(x,pad=pad_gene),axis=1)    
    return seqmap

def ariba_parse_stderr(stderr_fn):
    msgs = []
    summary = None
    with open(stderr_fn,'r') as inp:
        for line in inp:
            m = re.match(r'(WARNING:|ERROR:)(?:\sIn cluster log\s+([^:]+):)?\s*(.*)',
                line.strip())
            if m:
                tag,cluster_log,msg = m.groups()
                if cluster_log is None:
                    cluster_log = "main"
                msgs.append(dict(tag=tag,
                    cluster_log=cluster_log,
                    msg=msg))
    if len(msgs) > 0:
        msgs = pd.DataFrame(msgs)
        msgs_gr = msgs.groupby(["tag","msg"])["cluster_log"].\
            apply(lambda x: ' in ' + ','.join(x)).\
            reset_index(level=['tag', 'msg'])
        ## consider earlier messages as more informative
        msg_order = msgs[['tag','msg']].drop_duplicates().copy()
        msg_order["i_msg"] = np.arange(msg_order.shape[0])
        msg_gr = pd.merge(msg_order,msgs_gr)
        msg_gr = msg_gr.sort_values("i_msg")
        del msg_gr["i_msg"]
        msg_paste = lambda x: "{} {} {}".format(x[0],x[1],x[2])
        summary = msgs_gr.apply(msg_paste, axis=1, reduce=True).str.cat(sep=";")
    return summary


def extract_contigs(sample_id,report,assembled_genes,assembled_seqs,assemblies,
                    seq_map,ariba_stdout,ariba_stderr,
                    basecov_asm,basecov_ref,
                    report_out,fasta_out,status_out,
                    cut_to_ref=False,
                    pad_assembled=200,
                    pad_gene=200,
                    one_seq_per_contig=False,
                    sig_inp=None):
    """
    Extract contigs and aggregated manifest from Ariba outputs
    :param sample_id:
    :param report:
    :param assembled_genes:
    :param assembled_seqs:
    :param assemblies:
    :param report_out:
    :param fasta_out:
    :param cut_to_ref: Try to map assembled_seq to assemblies (should be exact substring) and then
    cut and pad by pad_assembled - full assemblies might be too long for RADAR.
    If not found, just use the full assembly. This is set to False by default because novel
    sequences can be trimmed by this operation.
    :param pad_assembled:
    :param sig_inp: Optional comma-separated string of input file cryptographic signatures. Each will result in a field `Asm_sig_inp1`, 'Asm_sig_inp2' etc
    """
    #
    aux_columns = \
    ['ContigSeqID',
     'ref_name',
     'gene',
     'var_only',
     'flag',
     'reads',
     'cluster',
     'ref_len',
     'ref_base_assembled',
     'pc_ident',
     'ctg',
     'ctg_len',
     'ctg_cov',
     'seq_len',
     '#ariba_ref_name',
     'SeqStatus']
    main_columns = \
    [
        'SampleID',
        'SeqID',
    ]
    ## prepare data for the input signature fields
    if sig_inp:
        sig_inp = OrderedDict([ ("Asm_sig_inp{}".format(i_sig+1),sig.strip()) for (i_sig,sig) in enumerate(sig_inp.strip().split(",")) ])
    else:
        sig_inp = OrderedDict()
    #aux_columns = aux_columns + list(sig_inp.keys())
    ## only leave those columns that are properly aggregated
    aggr_col_names = main_columns + aux_columns

    seqmap = load_seqmap_with_seqs(seq_map=seq_map,
        assembled_genes=assembled_genes,
        assembled_seqs=assembled_seqs,
        assemblies=assemblies)
    rep = pd.read_table(report)    
    if rep.shape[0] > 0:
        rep["SampleID"] = sample_id
        ## Replace dash. RADAR does not tolerate dash in seqid because it uses it to separate the added suffix in FASTA files
        seq_id_root = sample_id.replace("-","_")
        ## Having some better better aggr for the flag would be useful
        ## Note that "gene=1" can get replaced when multirecords are present (samp 2090)
        rep_aggr = rep.groupby(["ref_name","ctg"],as_index=False).head(1).copy()
        rep_aggr["rep_iid"] = np.arange(rep_aggr.shape[0])
        ## Load base coverage and overwrite the report field ctg_cov with the median base coverage
        basecov = pd.read_table(basecov_asm)
        basecov.rename({"#RefName":"ctg","Coverage":"ctg_cov"},axis="columns",inplace=True)
        ctg_cov = basecov.groupby("ctg",as_index=False).ctg_cov.median()
        del rep_aggr["ctg_cov"]
        rep_aggr = pd.merge(rep_aggr,ctg_cov,on="ctg",how="left")
        assert not rep_aggr.ctg_cov.isnull().any(), "Some contigs are absent from assembly base coverage table"
        if one_seq_per_contig:
            seqmap = prune_seqmap(seqmap)
        ##TODO: we need to align each match to the reference with nucmer or minimap2 here and flip
        ## if necessary. It is possible that the match is oriented in revcomp direction to the
        ## reference if several shorter matches on the contig influenced the orientation together.
        seqmap = seqmap_extend_matches(seqmap,pad_assembled=pad_assembled,pad_gene=pad_gene)
        ## get a single gene records per ref,contig,hit_id for later merging. should not be necessary
        seqmap_uniq_gene = seqmap.loc[seqmap.seq_type=="gene",["ref_name","contig","seq_type","hit_id"]].drop_duplicates()
        ## left-merge with seqmap to get a matching gene (if any) for each match
        seqmap = pd.merge(seqmap,seqmap_uniq_gene,on=["ref_name","contig","hit_id"],how="left",suffixes=["","_gene"])
        ## seed final padded sequence table with all genes
        rep_pad = pd.merge(rep_aggr,seqmap.loc[seqmap.seq_type=="gene",["ref_name","contig","seq_id","seq_pad","seq_type"]],
            left_on=["ref_name","ctg"],right_on=["ref_name","contig"])
        if cut_to_ref:
            ## append to rep_pad all matches that do not have a matching gene
            rep_match = pd.merge(rep_aggr,seqmap.loc[np.logical_and(seqmap.seq_type=="match",seqmap.seq_type_gene.isnull()),
                ["ref_name","contig","seq_id","seq_pad","seq_type"]],
                left_on=["ref_name","ctg"],right_on=["ref_name","contig"])
            rep_pad = pd.concat([rep_pad,rep_match],axis=0)
        ## append to rep_pad all contigs if there is not record yet in rep_pad from that contig
        rep_ctg = pd.merge(rep_aggr.loc[-rep_aggr.rep_iid.isin(rep_pad.rep_iid)],
            seqmap.loc[seqmap.seq_type=="asm",["contig","seq_id","seq_pad","seq_type"]],
            left_on=["ctg"],right_on=["contig"],how="left")
        rep_pad = pd.concat([rep_pad,rep_ctg],axis=0)
        ## rename fields for final output and generate serial SeqID
        rep_pad.rename(dict(seq_type="SeqStatus",seq_id="ContigSeqID"),axis="columns",inplace=True)
        rep_pad.loc[rep_pad.SeqStatus.isnull(),"SeqStatus"] = "miss"
        rep_pad["rep_iid"] = np.arange(1,rep_pad.shape[0]+1)
        rep_pad["SeqID"] = rep_pad.apply(lambda rec,seq_id_root=seq_id_root: "{}_{}".format(seq_id_root,rec.rep_iid),axis=1)
        rep_pad["seq_len"] = rep_pad.seq_pad.str.len()
    else:
        rep_pad = pd.DataFrame(columns=aggr_col_names)

    with open(fasta_out,"w") as fa_out:
        for r_rep in rep_pad.itertuples():
            if r_rep.SeqStatus != "miss":
                rec_out = SeqRecord(id=r_rep.SeqID,seq=Seq(r_rep.seq_pad,ambiguous_dna),description="")
                SeqIO.write([rec_out],fa_out,"fasta")
    
    rep_pad = rep_pad[aggr_col_names]
    renames = { nm : "Asm_" + nm.lstrip("#") for nm in aux_columns }
    del renames["flag"] # ariba expandflag will need its original name
    rep_pad.rename(columns=renames,inplace=True)
    rep_pad.to_csv(report_out, index=False, sep="\t")
    ariba_warn = ariba_parse_stderr(ariba_stderr)
    if not ariba_warn:
        ariba_warn = "Finished"
    df_status = pd.DataFrame([dict(SampleID=sample_id,
        Asm_Msg=ariba_warn[:128].replace("\t"," "))])
    for sig_fld in sig_inp:
        df_status[sig_fld] = sig_inp[sig_fld]

    df_status.to_csv(status_out,index=False, sep="\t")

def combine_samples(sample_results,manifest,manifest_out,fasta_out):
    import pandas as pd
    from Bio import SeqIO
    import subprocess
    man = workflow_util.load_manifest(manifest)
    assert not man.duplicated("SampleID").any(), "SampleID in the manifest contains duplicated values"
    sres = pd.read_table(sample_results)
    rep_all = []
    status_all = []
    with util.open_text_py23(fasta_out,"w") as fout:
        for r_sres in sres.itertuples():
            rep = workflow_util.load_manifest(r_sres.report)
            rep_all.append(rep)
            status = workflow_util.load_manifest(r_sres.status)
            status_all.append(status)
            with util.open_text_py23(r_sres.sequence,"r") as finp:
                SeqIO.write(SeqIO.parse(finp,"fasta"),fout,"fasta")
    rep_all = pd.concat(rep_all, axis=0).reset_index(drop=True)
    status_all = pd.concat(status_all, axis=0).reset_index(drop=True)    
    if rep_all.shape[0] > 0:
        # spaces here break output of `ariba expandflag`
        rep_all["free_text"] = "."
        rep_all_f = util.make_work_file(manifest_out,location="cwd")
        rep_all_f_exp = rep_all_f+".exp"
        rep_all.to_csv(rep_all_f, index=False, sep="\t")
        subprocess.check_call(["ariba","expandflag",rep_all_f,rep_all_f_exp])
        rep_all = workflow_util.load_manifest(rep_all_f_exp)
    rep_all.rename(columns=dict(flag="Asm_flag"),inplace=True)
    if "free_text" in list(rep_all):
        rep_all.drop(columns=["free_text"],inplace=True)
    rep_not_man = set(rep_all.SampleID) - set(man.SampleID)
    if len(rep_not_man) > 0:
        raise ValueError("Some SampleID values in the reports are not in the manifest: {}".format(rep_not_man))
    rep_mrg = pd.merge(man,rep_all,on=["SampleID"],how="left",suffixes=["_Manifest",""])
    rep_mrg.reset_index(drop=True,inplace=True)
    ## this will never happen in the left merge, but won't hurt to check
    assert set(man.SampleID) == set(rep_mrg.SampleID), "Some SampleID values in manifest are lost after merging with the report"
    assert num_util.series_equal(man.SampleID,status_all.SampleID), "Expected 1-to-1 relations on SampleID between manifest and status files"    
    rep_mrg = pd.merge(rep_mrg,status_all,on=["SampleID"],suffixes=["_Manifest",""])
    rep_mrg.reset_index(drop=True,inplace=True)    
    assert num_util.series_equal(man.SampleID,rep_mrg.SampleID.drop_duplicates()), "Expected 1-to-many relations on SampleID between input and output manifest files"
    rep_mrg.to_csv(manifest_out, index=False, sep="\t", na_rep="NA")


def _cwl_file(path,to_abspath=True):
    log.debug(path)
    if to_abspath:
        path = os.path.abspath(path)
    return { 'class' : 'File', 'path' : path}

@argh.arg("--cwl-inputs",help="The initial workflow inputs file. Missing keys will be populated from other arguments or defaults")
@argh.arg("--datadir",help="Optional top directory of the sequence files in the manifest. "
                           "If missing, manifest should list absolute paths")
@argh.arg("--allowed-roots",help="Optional comma-separated list of directories that are allowed roots of the sequence files in the manifest. "
                           "If provided, each manifest file should be located under one of those.")
@argh.arg("--micgent-data",help="Directory with MICGENT database files. Required if adapter and spikes files are not present in "
    "the --cwl-inputs")
@argh.arg("--assembly-policy",help="Sets default values for several parameters (unless those are defined directly) "
                                   " in order to meet some common use cases. Any affected parameters that are already "
                                   " defined in the input config will be kept as-is.",
          choices=["wgs_fermilite","wgs_spades","amplicon"])
@argh.arg("--only-configure-inputs",help="Fill in the final workflow inputs file and exit")
## need type=int because None is a default and argh cannot deduce the correct type from function signature
@argh.arg("--sor-sfx-version",help="SeqID FASTA version for system of record (SOR)",type=int)
@argh.arg("--sig-inp-key",help="Key to use when computing signatures of the input files. It has to be kept constant across re-runs for signature comparisons.")
def run_extraction_wf(cwl_inputs=None,
                      datadir="",
                      allowed_roots=None,
                      micgent_data=None,
                      assembly_policy=None,
                      prepareref_tgz=None,
                      manifest=None,
                      primer_literals=None,
                      ref_common=None,
                      micgentjs_tgz=None,
                      filter_asm_args=None,
                      filter_asm_default_ctg_len_min=0,
                      filter_asm_default_ctg_cov_min=20,
                      filter_asm_default_ctg_cov_ratio_min=0.05,
                      deterministic=False,
                      sor_pack_out=None,
                      seq_out=None,
                      manifest_out=None,
                      web_tar_out=None,
                      web_dir_out=None,
                      sor_sfx_target=None,
                      sor_sfx_version=None,
                      sig_inp_key=None,                      
                      debug=False,
                      cwl_runner_name="toil",
                      cwl_runner_config=None,
                      outdir=None,
                      only_configure_inputs=False):
    """
    Prepare and execute gene extractor workflow.

    This function will finalize the provided, potentially incomplete, YAML workflow input file
    using sensible defaults and built-in database files when necessary, execute the CWL workflow
    with the requested runner, and then package and copy output files to the requested locations.

    :param prepareref_tgz:
    :param manifest:
    :param primer_literals:
    :param ref_common:
    :param micgentjs_tgz:
    :param filter_asm_args:
    :param filter_asm_default_ctg_len_min:
    :param sor_pack_out:
    :param seq_out:
    :param manifest_out:
    :param web_tar_out:
    :param web_dir_out:
    :param sor_sfx_target:
    :param sor_sfx_version:
    :param cwl_runner_name:
    :param cwl_runner_config:
    :param outdir:
    """
    pjoin = os.path.join
    pkg_data = resources.get_pkg_data_dir("MICGENT")
    cwl_dir = os.path.join(pkg_data,"cwl")
    assert cwl_runner_name in ("toil",),"Only Toil CWL runner is currently supported"
    out_inputs = yaml_util.load_yaml(cwl_inputs) if cwl_inputs else {}
    curdir = os.getcwd()
    if micgent_data:
        micgent_data = os.path.abspath(micgent_data)
    if not outdir:
        outdir = os.path.join(curdir,"out")
    else:
        outdir = os.path.abspath(outdir)
    if not os.path.exists(outdir):
        os.makedirs(outdir)
    if web_dir_out:
        web_dir_out = os.path.abspath(web_dir_out)
    def _define_out_file(inp,def_base,def_dir):
        if not inp:
            inp = os.path.join(def_dir,def_base)
        return os.path.abspath(inp)
    def _db_file(rel_path):
        assert micgent_data,"Path to MICGENT DB must be defined"
        return _cwl_file(pjoin(micgent_data, rel_path))
    def _setdefault_db_file(key,rel_path):
        if key not in out_inputs:
            out_inputs[key] = _db_file(rel_path)
        return out_inputs[key]['path']
    def _set_val(key,val,default=None,converter=lambda x:x):
        if val is None:
            val = out_inputs.get(key,None)
        if val is None:
            val = default
        if val is not None:
            out_inputs[key] = converter(val)
        return out_inputs.get(key,None)
    def _set_file(key, path, default=None):
        ret = _set_val(key,path,default,converter=_cwl_file)
        if ret is not None:
            return ret["path"]
        return None

    sor_pack_out = _define_out_file(sor_pack_out,"sor_pack.tgz",outdir)
    seq_out = _define_out_file(seq_out, "seq_out.fasta", outdir)
    manifest_out = _define_out_file(manifest_out, "manifest_out.tsv", outdir)
    web_tar_out = _define_out_file(web_tar_out, "web.tar", outdir)

    manifest = _set_file('manifest',manifest)

    samples = os.path.join(outdir,"gene_extractor_samples.yaml")
    workflow_util.manifest_to_cwl(manifest,samples,datadir=datadir,allowed_roots=allowed_roots)

    prepareref_tgz = _set_file('prepareref_tgz',prepareref_tgz)

    if primer_literals:
        out_inputs['primer_literals'] = primer_literals.split(",")

    ref_common = _set_file('ref_common',ref_common)

    micgentjs_tgz = _set_file('micgentjs_tgz',micgentjs_tgz)

    if not micgentjs_tgz:
        micgentjs_src_dir = os.path.join(pkg_data,"web","micgentjs")
        micgentjs_tgz = os.path.join(curdir,"micgentjs.tgz")
        util.dir_to_tar(micgentjs_src_dir,micgentjs_tgz,with_dir=True)
        micgentjs_tgz = _set_file('micgentjs_tgz', micgentjs_tgz)

    _setdefault_db_file("spikes_file","qa/bbtools/phiX.fa")
    _setdefault_db_file("adapter_file", "qa/bbtools/bbduk_adapters.fa")

    _set_val('filter_asm_args',filter_asm_args,
             '{policy: default, ctg_len_min: %s, ctg_cov_ratio_min: %s, ctg_cov_min: %s}' % \
             (filter_asm_default_ctg_len_min,filter_asm_default_ctg_cov_ratio_min,filter_asm_default_ctg_cov_min))

    _set_val('sor_sfx_target',sor_sfx_target)

    _set_val('sor_sfx_version',sor_sfx_version)

    _set_val('sig_inp_key',sig_inp_key)
    
    assembly_policy = _set_val('assembly_policy',assembly_policy)

    if assembly_policy:
        if assembly_policy == "amplicon":
            out_inputs.setdefault("assembly_cov",100000000)
            out_inputs.setdefault("assembler","plugin")
            out_inputs.setdefault("plugin_asm_options","python -m MICGENT.ariba_asm_plugin asm-for-skewed-coverage")
        elif assembly_policy == "wgs_spades":
            out_inputs.setdefault("assembler","spades")
            out_inputs.setdefault("assembly_cov",100)
        elif assembly_policy == "wgs_fermilite":
            out_inputs.setdefault("assembler", "fermilite")
            out_inputs.setdefault("assembly_cov", 100)
        else:
            raise ValueError("Unknown assembly policy requested: {}".format(assembly_policy))

    if not deterministic:
        deterministic = None # so that default pulling logic works
    _set_val('deterministic',deterministic,False)

    if not debug:
        debug = None # so that default pulling logic works
    _set_val('debug',debug,False)

    out_inputs['samples'] = {}
    out_inputs['samples']['$import'] = samples

    wf_inputs = os.path.join(outdir,"gene_extractor.yaml")

    ## pass workflow inputs as another file input for provenance keeping
    out_inputs['wf_inputs'] = _cwl_file(wf_inputs)

    yaml_util.dump_yaml(out_inputs,wf_inputs)

    if not only_configure_inputs:

        wf = os.path.join(cwl_dir, "gene_extractor.cwl")

        runner_func_name = "run_{}".format(cwl_runner_name)
        runner_func = getattr(cwl_runner,runner_func_name)
        ## cwl_runner_config is expected to be produced by arg_parsing.dump_sig, and be keyed by the function name
        runner_conf = yaml_util.load_yaml(cwl_runner_config)[runner_func_name] if cwl_runner_config else {}
        runner_func(wf=wf,wf_inputs=wf_inputs,**runner_conf)

        def _copy_if_diff_name(finp,fout):
            if not finp is None:
                if not (os.path.exists(fout) and os.path.samefile(finp,fout)):
                    shutil.copy(finp,fout)
        with util.chdir(outdir):
            _copy_if_diff_name("sor_pack.tgz", sor_pack_out)
            _copy_if_diff_name("seq_out.fasta", seq_out)
            _copy_if_diff_name("manifest_out.tsv", manifest_out)
            _copy_if_diff_name("web.tar", web_tar_out)
            _copy_if_diff_name(prepareref_tgz,"prepareref.tgz")
            if web_dir_out:
                ##that checks that dir either does not exits or empty, and creates intermediate dirs
                util.tar_to_dir("web.tar",web_dir_out)


def _validate_radar_timestamp(timestamp):
    if not (len(timestamp) == 12 and timestamp.isdigit()): 
        log.warning("RADAR timestamp does not have a legitimate format: {}".format(timestamp))
        return False
    try:
        time.strptime('201810161812','%Y%m%d%H%M')
    except ValueError:
        log.exception("RADAR timestamp failed validation parsing")
        return False
    return True

@argh.arg("--jar",help="Java JAR file that implements RADAR upload CLI")
@argh.arg("--rmi-registry",help="Java RMI URL to use as endpoint for RADAR upload")
def radar_upload(res_pack,stage,department,project_code,study_code,timestamp,user_email,
                galaxy_provenance_dir,
                jar=None,
                rmi_registry=None):
    assert jar and rmi_registry, "jar and rmi_registry arguments must be provided"
    stage = stage.strip()
    department = department.strip()
    project_code = project_code.strip()
    study_code = study_code.strip()
    timestamp = timestamp.strip()
    user_email = user_email.strip()
    assert stage and department and project_code and study_code and timestamp, \
        "Need all non-empty components for the upload RADAR manifest file name"
    assert user_email, "Need a defined user for RADAR upload: {}".format(user_email)
    assert _validate_radar_timestamp(timestamp), "Timestamp failed format validation: {}".format(timestamp)
    jar = os.path.abspath(jar)
    assert os.path.isdir(galaxy_provenance_dir), "Galaxy provenance parameter must be an existing directory"
    galaxy_provenance_dir = os.path.abspath(galaxy_provenance_dir)
    pack_dir = "res_pack"
    util.tar_to_dir(res_pack,pack_dir)
    user = user_email.split("@")[0].strip()
    assert user, "Need a defined user for RADAR upload: {}".format(user)
    with util.chdir(pack_dir):
        radar_manifest = "./{}_{}_{}_{}_{}.txt".format(stage,department,project_code,study_code,timestamp)
        shutil.move("manifest.tsv",radar_manifest)
        provenance = "provenance.tgz"
        provenance_dir = "provenance"
        util.tar_to_dir(provenance,provenance_dir)
        shutil.copytree(galaxy_provenance_dir,os.path.join(provenance_dir,"galaxy_provenance"))
        os.remove(provenance)
        toc_fn = os.path.join(provenance_dir,"toc.tsv")
        toc = pd.read_csv(toc_fn,dialect="excel-tab",dtype=dict(Index=str))
        toc = toc.append([{"Index":"9.1","Description":"Galaxy provenance data","File":"galaxy_provenance/job_user.yaml"}],
            ignore_index=True)
        toc.to_csv(toc_fn, index=False, sep="\t")        
        util.dir_to_tar(provenance_dir,provenance,with_dir=True)
        try:
            subprocess.check_call(["java",
                                   "-jar",
                                   jar,
                                   rmi_registry,
                                   user,
                                   radar_manifest,
                                   provenance])
        ## we need to supress all exception traces in check_call to make sure the RADAR RMI is not revealed in stderr
        except:
            exc_type, exc_value, exc_traceback = sys.exc_info()
            sys.exit("Call to RADAR upload interface has failed with the exception of type: {}".format(getattr(exc_type,'__name__',repr(exc_type))))

@argh.arg("--jar",help="Java JAR file that implements RADAR upload CLI")
@argh.arg("--rmi-registry",help="Java RMI URL to use as endpoint for RADAR upload")
def radar_check_user(user,
                 jar=None,
                 rmi_registry=None):
    assert jar and rmi_registry, "jar and rmi_registry arguments must be provided"
    jar = os.path.abspath(jar)
    user = user.strip()
    assert user, "Need a defined user for RADAR upload: {}".format(user)
    subprocess.check_call(["java",
                       "-jar",
                       jar,
                       "-v",
                       rmi_registry,
                       user],
                       stderr=subprocess.STDOUT)

## import package module and add argh entry points

def _main():
    parser = arg_parsing.ArghParserChainedConfig()
    parser.add_commands([
        extract_contigs,
        combine_samples,
        run_extraction_wf,
        radar_upload,
        radar_check_user
    ])
    parser.dispatch()


if __name__ == "__main__":
    _main()
