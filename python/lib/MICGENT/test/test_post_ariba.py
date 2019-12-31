from __future__ import unicode_literals
from __future__ import print_function
from __future__ import division
from __future__ import absolute_import
from future import standard_library
standard_library.install_aliases()
from builtins import *

from MICGENT import yaml_util
from MICGENT import arg_parsing
from MICGENT import cwl_runner
from MICGENT import gene_extractor_bench as ge_bench
from MICGENT import resources
from MICGENT import seq_util

import helpers

import pytest
from subprocess import check_call, check_output
import os
import shlex
import shutil

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet.IUPAC import ambiguous_dna
from Bio.Seq import Seq

import pandas as pd
import numpy as np

from MGT.Logging import *

pytestmark = pytest.mark.usefixtures("goto_cleandir_test","get_test_data_dir")

test_data = "test_data"

pjoin = os.path.join


def _compare_tables(expected_json,returned_manifest,pars,descr):
    """Load a tab-delimited table file and compare with the expected DataFrame passed as df.to_json(ordering='columns') results.
    Only columns from the expected will be considered, but the datatypes will be forced to the types from the returned."""
    man = pd.read_table(returned_manifest)
    expected_df = pd.read_json(expected_json,dtype=dict(man.dtypes))
    returned_df = man[expected_df.columns]
    assert expected_df.equals(returned_df), \
        "Test objective: {}. Unexpected values of the output manifest {} returned after testing with parameters: {}. Expected: {}. Returned: {}".\
        format(descr,returned_manifest,pars,expected_df.to_json(orient="columns"),returned_df.to_json(orient="columns"))


@pytest.mark.parametrize("descr,filter_yaml,manifests", [
            (
            "All pass",
            "{policy: default, ctg_len_min: 100, "
            "ctg_cov_ratio_min: 0.05, ctg_cov_min: 20, "
            "ctg_ref_bases_min: 100, ctg_ref_bases_ratio_min: 0.1}",
            [
                ("manifest_out.tsv",'{"Asm_SeqCountPostFilter":{"0":2,"1":2},"SampleID":{"0":"SASMA1927338_1","1":"SASMA1927338_1"},"SeqID":{"0":"SASMA1927338_1_1","1":"SASMA1927338_1_2"}}'),
                ("manifest_out_all.tsv",'{"Asm_Dropped":{"0":false,"1":false},"Asm_DroppedMultiseq":{"0":false,"1":false},"Asm_SeqCountPostFilter":{"0":2,"1":2},\
                    "SampleID":{"0":"SASMA1927338_1","1":"SASMA1927338_1"},"SeqID":{"0":"SASMA1927338_1_1","1":"SASMA1927338_1_2"}}'),
                ("manifest_out_sum.tsv",'{"Asm_Failed":{"0":false},"Asm_SeqCountFinal":{"0":2},"Asm_DroppedCount":{"0":0},\
                    "Asm_DroppedMultiseq":{"0":false},"Asm_SeqCountPostFilter":{"0":2},"SampleID":{"0":"SASMA1927338_1"}}')
            ]
            ),

            (
            "Multiseq Rejects",
            "{policy: default, ctg_len_min: 100, "
            "ctg_cov_ratio_min: 0.05, ctg_cov_min: 20, "
            "ctg_ref_bases_min: 100, ctg_ref_bases_ratio_min: 0.1, multiseq_policy: reject}",
            [
                ("manifest_out.tsv",'{"SampleID":{},"SeqID":{},"Asm_SeqCountPostFilter":{}}'),
                ("manifest_out_all.tsv",'{"Asm_Dropped":{"0":true,"1":true},"Asm_DroppedMultiseq":{"0":true,"1":true},"Asm_SeqCountPostFilter":{"0":2,"1":2},\
                    "SampleID":{"0":"SASMA1927338_1","1":"SASMA1927338_1"},"SeqID":{"0":"SASMA1927338_1_1","1":"SASMA1927338_1_2"}}'),
                ("manifest_out_sum.tsv",'{"Asm_Failed":{"0":true},"Asm_SeqCountFinal":{"0":0},"Asm_DroppedCount":{"0":2},\
                    "Asm_DroppedMultiseq":{"0":true},"Asm_SeqCountPostFilter":{"0":2},"SampleID":{"0":"SASMA1927338_1"}}')
            ]
            ),

            (
            "One dropped by contig length",
            "{policy: default, ctg_len_min: 790, "
            "ctg_cov_ratio_min: 0.05, ctg_cov_min: 20, "
            "ctg_ref_bases_min: 100, ctg_ref_bases_ratio_min: 0.1}",
            [
                ("manifest_out.tsv",'{"Asm_SeqCountPostFilter":{"0":1},"SampleID":{"0":"SASMA1927338_1"},"SeqID":{"0":"SASMA1927338_1_1"}}'),
                ("manifest_out_all.tsv",'{"Asm_Dropped":{"0":false,"1":true},"Asm_DroppedMultiseq":{"0":false,"1":false},"Asm_SeqCountPostFilter":{"0":1,"1":1},\
                    "SampleID":{"0":"SASMA1927338_1","1":"SASMA1927338_1"},"SeqID":{"0":"SASMA1927338_1_1","1":"SASMA1927338_1_2"}}'),
                ("manifest_out_sum.tsv",'{"Asm_Failed":{"0":false},"Asm_SeqCountFinal":{"0":1},"Asm_DroppedCount":{"0":1},\
                    "Asm_DroppedMultiseq":{"0":false},"Asm_SeqCountPostFilter":{"0":1},"SampleID":{"0":"SASMA1927338_1"}}')
            ]
            ),

            (
            "One dropped by ref coverage ratio",
            "{policy: default, ctg_len_min: 100, "
            "ctg_cov_ratio_min: 0.05, ctg_cov_min: 20, "
            "ctg_ref_bases_min: 100, ctg_ref_bases_ratio_min: 0.8}",
            [
                ("manifest_out.tsv",'{"Asm_SeqCountPostFilter":{"0":1},"SampleID":{"0":"SASMA1927338_1"},"SeqID":{"0":"SASMA1927338_1_1"}}'),
                ("manifest_out_all.tsv",'{"Asm_Dropped":{"0":false,"1":true},"Asm_DroppedMultiseq":{"0":false,"1":false},"Asm_SeqCountPostFilter":{"0":1,"1":1},\
                    "SampleID":{"0":"SASMA1927338_1","1":"SASMA1927338_1"},"SeqID":{"0":"SASMA1927338_1_1","1":"SASMA1927338_1_2"}}'),
                ("manifest_out_sum.tsv",'{"Asm_Failed":{"0":false},"Asm_SeqCountFinal":{"0":1},"Asm_DroppedCount":{"0":1},\
                    "Asm_DroppedMultiseq":{"0":false},"Asm_SeqCountPostFilter":{"0":1},"SampleID":{"0":"SASMA1927338_1"}}')
            ]
            ),

            (
            "One dropped by contig length, ref coverage ratio and length",
            "{policy: default, ctg_len_min: 800, "
            "ctg_cov_ratio_min: 0.05, ctg_cov_min: 20, "
            "ctg_ref_bases_min: 900, ctg_ref_bases_ratio_min: 0.8}",
            [
                ("manifest_out.tsv",'{"SampleID":{},"SeqID":{},"Asm_SeqCountPostFilter":{}}'),
                ("manifest_out_all.tsv",'{"Asm_DroppedMultiseq":{"0":false,"1":false},"Asm_SeqCountPostFilter":{"0":0,"1":0},\
                    "SampleID":{"0":"SASMA1927338_1","1":"SASMA1927338_1"},"SeqID":{"0":"SASMA1927338_1_1","1":"SASMA1927338_1_2"}}'),
                ("manifest_out_sum.tsv",'{"Asm_Failed":{"0":true},"Asm_SeqCountFinal":{"0":0},"Asm_DroppedCount":{"0":2},\
                    "Asm_DroppedMultiseq":{"0":false},"Asm_SeqCountPostFilter":{"0":0},"SampleID":{"0":"SASMA1927338_1"}}')
            ]
            )
            ]
    )
def test_filter_assemblies(descr,filter_yaml,manifests):
    test_data = os.path.abspath(globals()["test_data"])
    filter_inp_dir = pjoin(test_data,"gene_extractor/post_ariba/filter")
    with helpers.mkchdir("filter_assemblies"):
        cmd = ("python -m MICGENT.post_ariba filter-assemblies "
            "--args '{filter_yaml}' "
            "--manifest {filter_inp_dir}/manifest_out.tsv "
            "--contigs {filter_inp_dir}/seq_out.fasta").format(**locals())
        print(cmd)
        check_call(shlex.split(cmd))
        for returned_manifest,expected_json in manifests:
            _compare_tables(expected_json,returned_manifest,filter_yaml,descr)

@pytest.mark.skip(reason="Need to provide basecov-asm and basecov-ref positional arguments")
def test_extract_contigs(request):
    test_data = os.path.abspath(globals()["test_data"])
    inp_dir = pjoin(test_data,"gene_extractor/extract_contigs")
    with helpers.mkchdir("extract_contigs"):
        cmd = ("python -m MICGENT.gene_extractor extract-contigs "
            "--cut-to-ref --pad-assembled 200 --sig-inp 123,123 "
            "SMA1828869 "
            "{inp_dir}/report.tsv "
            "{inp_dir}/assembled_genes.fa.gz "
            "{inp_dir}/assembled_seqs.fa.gz "
            "{inp_dir}/assemblies.fa.gz "
            "{inp_dir}/seq_map.tsv "
            "{inp_dir}/stdout.log "
            "{inp_dir}/stderr.log "
            "report_out.tsv seq_out.fasta status.tsv").format(**locals())
        print(cmd)
        check_call(shlex.split(cmd))
        man = pd.read_table("report_out.tsv")
        assert man.shape[0] == 2
        assert man[man.SeqStatus=="gene"].shape[0] == 1
        assert man[man.SeqStatus=="match"].shape[0] == 1

@pytest.mark.skip(reason="Need to provide basecov-asm and basecov-ref positional arguments")
def test_extract_contigs_stop_codon_revcomp(request):
    """Gene found with a stop codon and overlaps the assembled match partially,
    and found either on positive or on negative strand.
    Test that we correctly identify the correspodence between the gene and the match,
    and will not output the match, and also output the gene on the correct strand"""
    test_data = os.path.abspath(globals()["test_data"])
    for extr_cont_dir in ("extract_contigs_stop_codon","extract_contigs_stop_codon_revcomp"):
        inp_dir = pjoin(test_data,"gene_extractor",extr_cont_dir)
        with helpers.mkchdir(extr_cont_dir):
            for (pad_assembled, pad_gene) in [(0,0),(1,1),(200,200),(200,0),(0,200)]:
                cmd = ("python -m MICGENT.gene_extractor extract-contigs "
                    "--cut-to-ref --pad-assembled {pad_assembled} --pad-gene {pad_gene} "
                    "SMA807030 "
                    "{inp_dir}/report.tsv "
                    "{inp_dir}/assembled_genes.fa.gz "
                    "{inp_dir}/assembled_seqs.fa.gz "
                    "{inp_dir}/assemblies.fa.gz "
                    "{inp_dir}/seq_map.tsv "
                    "{inp_dir}/stdout.log "
                    "{inp_dir}/stderr.log "
                    "report_out.tsv seq_out.fasta status.tsv").format(**locals())
                print(cmd)
                check_call(shlex.split(cmd))
                man = pd.read_table("report_out.tsv")
                assert man.shape[0] == 1
                assert man.loc[0,"SeqStatus"] == "gene"
                seq_out = seq_util.fasta_to_df("seq_out.fasta",seq_format_out="str")
                seq_gene = seq_util.fasta_to_df("{inp_dir}/assembled_genes.fa.gz".format(inp_dir=inp_dir),seq_format_out="str")
                assert seq_out.shape[0] == 1
                assert str(seq_gene.seq[0]) in str(seq_out.seq[0])

def test_extract_contigs_skewed_cov_filter_rsv(request):
    test_data = os.path.abspath(globals()["test_data"])
    inp_dir = pjoin(test_data,"gene_extractor/extract_contigs_skewed_cov_filter_rsv")
    with helpers.mkchdir("extract_contigs"):
        cmd = ("python -m MICGENT.gene_extractor extract-contigs "
            "--cut-to-ref --pad-assembled 200 --sig-inp 123,123 "
            "SMA589764 "
            "{inp_dir}/report.tsv "
            "{inp_dir}/assembled_genes.fa.gz "
            "{inp_dir}/assembled_seqs.fa.gz "
            "{inp_dir}/assemblies.fa.gz "
            "{inp_dir}/seq_map.tsv "
            "{inp_dir}/stdout.log "
            "{inp_dir}/stderr.log "
            "{inp_dir}/basecov_asm.txt "
            "{inp_dir}/basecov_ref.txt "
            "report_out.tsv seq_out.fasta status.tsv").format(**locals())
        print(cmd)
        check_call(shlex.split(cmd))
        man = pd.read_table("report_out.tsv")
        assert man.shape[0] == 1
        assert (man.Asm_ctg_cov == 10).all()
