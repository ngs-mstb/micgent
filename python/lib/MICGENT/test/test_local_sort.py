from MICGENT import sig

import helpers

import pytest

import os
import glob
import shlex
from os.path import join as pjoin
from subprocess import check_call

from MGT.Logging import *

pytestmark = pytest.mark.usefixtures("goto_cleandir_test","get_test_data_dir")

test_data = "test_data"
test_data_local_sort = os.path.join(test_data,"local_sort")

def test_sort_by_ref_in_window():
    from MICGENT.local_sort import sort_by_ref_in_window
    inp_ref = [1,2,3,4,5,6,7,8,9]
    inp_que = [2,3,6,4,7,8,9]
    out_exp = [2,3,4,6,7,8,9]
    window_size = 2
    out_obs = []
    for el_id, el_body in sort_by_ref_in_window(inp_ref=inp_ref,
        inp_que=zip(inp_que,inp_que),
        window_size=window_size):
        out_obs.append(el_id)
        #print([el_id,el_body])
    assert(out_obs == out_exp)

def test_fastq_sort_by_ref_in_window():
    from MICGENT.local_sort_fastq import fastq_sort_by_ref_in_window
    out_fq_fn = "interleaved.sorted.fq"
    with open(pjoin(test_data_local_sort,"interleaved.fastq"),"r") as inp_que,\
        open(pjoin(test_data_local_sort,"interleaved.ref.id"),"r") as inp_ref,\
        open(out_fq_fn,"w") as out:
        fastq_sort_by_ref_in_window(inp_ref=inp_ref,inp_que=inp_que,window_size=4,
        out=out)
    with open(out_fq_fn,"r") as inp_out:
        out_body_str = inp_out.read()
    with open(pjoin(test_data_local_sort,"interleaved.sorted.fq"),"r") as inp_exp:
        exp_body_str = inp_exp.read()
    assert(out_body_str==exp_body_str)

def test_fastq_sort_by_ref_in_window_cli():
    inp_que_fn = pjoin(test_data_local_sort,"interleaved.fastq")
    inp_ref_fn = pjoin(test_data_local_sort,"interleaved.ref.id")
    out_fq_fn = "interleaved.sorted.fq"
    cmd = ("cat {inp_que_fn} | python -m MICGENT.local_sort_fastq fastq-sort-by-ref-in-window " 
        "--inp-ref {inp_ref_fn} "
        "--out {out_fq_fn} "
        "--window-size 4").format(**locals()) 
    #check_call(shlex.split(cmd))
    check_call(cmd,shell=True)
    with open(out_fq_fn,"r") as inp_out:
        out_body_str = inp_out.read()
    with open(pjoin(test_data_local_sort,"interleaved.sorted.fq"),"r") as inp_exp:
        exp_body_str = inp_exp.read()
    assert(out_body_str==exp_body_str)
