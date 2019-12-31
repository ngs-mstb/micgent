from __future__ import unicode_literals
from __future__ import print_function
from __future__ import division
from __future__ import absolute_import
from future import standard_library
standard_library.install_aliases()
from builtins import *

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

def test_run_step():
    from MICGENT.pysteps import run_step
    run_step("(printenv MICGENT_TEST_OK; printenv) > tmp.log",env=dict(MICGENT_TEST_OK="OK"),shell=True)
    with open("tmp.log","r") as inp:
        line = inp.readline().strip()
        assert line == "OK"

@helpers.skip_no_conda
def test_run_step_conda_root():
    from MICGENT.pysteps import run_step
    run_step("conda env --help | grep -i conda",use_conda=True)


def clean_reads(request,target="rsv",full_size=False):
    if target == "rsv":
        data_dir = pjoin(request.config.getoption('--large-test-data'),"rsv","SA2017" if full_size else "SA2017_sub")
        inp_reads = [ os.path.join(data_dir,"SMA1646701_S1_L001_R{}_001.fastq.gz".format(i_read)) for i_read in (1,2) ]
    elif target == "sa":
        data_dir = pjoin(request.config.getoption('--large-test-data'),"sa","sast2")
        inp_reads = [ os.path.join(data_dir,"SMA1828869_S1_L001_R{}_001.fastq.gz".format(i_read)) for i_read in (1,2) ]        
    inp_reads1 = inp_reads[0]
    inp_reads2 = inp_reads[1]
    with helpers.mkchdir("clean_reads_rsv"):
        out_reads1 = "cleaned.1.fq"
        out_reads2 = "cleaned.2.fq"
        out_qc_reads_cmd = (" --out-qc-before-reads qc.before.1.fq --out-qc-before-reads2 qc.before.2.fq "
            "--out-qc-after-reads qc.after.1.fq --out-qc-after-reads2 qc.after.2.fq ")
        for extra_steps in ("",
                out_qc_reads_cmd,
                "--primer-literals AGTGTTCAAYTTYGTWCCYTG,YTACCATTCAAGCAATGACCTC",
                "--clumpify --filter-spikes --primer-literals AGTGTTCAAYTTYGTWCCYTG,YTACCATTCAAGCAATGACCTC",
                out_qc_reads_cmd + \
                "--clumpify --filter-spikes --primer-literals AGTGTTCAAYTTYGTWCCYTG,YTACCATTCAAGCAATGACCTC"
                            ):
            cmd = ("python -m MICGENT.pysteps clean-reads --inp-reads {inp_reads1} --inp-reads2 {inp_reads2} "
                "--out-reads {out_reads1} --out-reads2 {out_reads2} "
                "--threads 4 --out-stats stats.log --deterministic " + extra_steps).format(**locals()) 
            check_call(shlex.split(cmd))
            s1 = sig.file_sig(out_reads1)
            s2 = sig.file_sig(out_reads2)
            ## save first replicate in case test fails and we need to debug
            os.rename(out_reads1,out_reads1+".r1")
            os.rename(out_reads2,out_reads2+".r1")
            check_call(shlex.split(cmd))
            assert os.path.getsize(out_reads1) > 0, "Output file {} has zero size".format(out_reads1)
            assert os.path.getsize(out_reads2) > 0, "Output file {} has zero size".format(out_reads2)
            assert sig.file_sig_cmp(s1,out_reads1), "Output files differ between repeated runs: {}".format(cmd)
            assert sig.file_sig_cmp(s2,out_reads2), "Output files differ between repeated runs: {}".format(cmd)
            if "--out-qc" in cmd:
                out_qc_reads_file_sizes = [ os.path.getsize(f) for f in glob.glob("qc.*.[12].fq") ]
                assert all([ _ > 0 for _ in out_qc_reads_file_sizes ]), "Some subsampled read file have zero size"
                assert len(out_qc_reads_file_sizes)== 4, "Some subsampled read files are missing"
            for fq in glob.glob("*.fq"):
                os.remove(fq)


@pytest.mark.slow
@pytest.mark.micgent_data
@pytest.mark.large_test_data
def test_clean_reads_rsv(request):
    clean_reads(request,target="rsv",full_size=False)

@pytest.mark.slow
@pytest.mark.micgent_data
@pytest.mark.large_test_data
def test_clean_reads_full_size_rsv(request):
    clean_reads(request,target="rsv",full_size=True)

@pytest.mark.slow
@pytest.mark.micgent_data
@pytest.mark.large_test_data
def test_clean_reads_full_size_sa(request):
    """This is a contrived test because we are trimming RSV primers on the Staph reads.
    We just need large enough input files in order to test the Bash pipe CPU utilization"""
    clean_reads(request,target="sa",full_size=True)
