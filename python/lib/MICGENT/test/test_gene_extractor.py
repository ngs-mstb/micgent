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

import helpers

import pytest
from subprocess import check_call, check_output
import os
import shlex
import shutil
import pandas as pd
import numpy as np

from MGT.Logging import *

pytestmark = pytest.mark.usefixtures("goto_cleandir_test","get_test_data_dir")

test_data = "test_data"

pjoin = os.path.join

@pytest.mark.micgent_data
@pytest.mark.large_test_data
def test_run_extractor_config_gen_rsv(request):
    test_data = os.path.abspath(globals()["test_data"])
    micgent_data = request.config.getoption('--micgent-data')
    conf_dir = pjoin(test_data,"gene_extractor/rsv")
    data_dir = pjoin(request.config.getoption('--large-test-data'),"rsv/SA2017")
    with helpers.mkchdir("extractor_config_gen_rsv"):
        ## Run w/o existing workflow inputs
        cmd_base = ("python -m MICGENT.gene_extractor run-extraction-wf --micgent-data {micgent_data} "
            "--datadir {data_dir} --prepareref-tgz {conf_dir}/ref_AB_single.trimmed.arbref.tgz "
            "--manifest {conf_dir}/sample_manifest_two.txt --ref-common {conf_dir}/ref_AB_single.fasta "
            "--only-configure-inputs").format(**locals())
        cmd = cmd_base
        check_call(shlex.split(cmd))
        assert os.path.exists("micgentjs.tgz")
        wf_inp1 = "out/gene_extractor.yaml"
        assert os.path.exists(wf_inp1)
        wf_inp1 = yaml_util.load_yaml(wf_inp1)
        print(wf_inp1)
        assert os.path.basename(wf_inp1["spikes_file"]["path"]) == "phiX.fa"
        assert "primer_literals" not in wf_inp1
        assert "assembler" not in wf_inp1
        assert yaml_util.get_arg_as_yaml(wf_inp1["filter_asm_args"])["ctg_len_min"] == 0
        assert os.path.basename(wf_inp1["ref_common"]["path"]) == "ref_AB_single.fasta"
        assert os.path.basename(wf_inp1["micgentjs_tgz"]["path"]) == "micgentjs.tgz"

        ## Run again with existing workflow inputs
        cmd = cmd_base + " --cwl-inputs {conf_dir}/cwl_inputs.yaml".format(**locals())
        check_call(shlex.split(cmd))
        wf_inp2 = "out/gene_extractor.yaml"
        assert os.path.exists(wf_inp2)
        wf_inp2 = yaml_util.load_yaml(wf_inp2)
        print(wf_inp2)
        assert wf_inp2["primer_literals"][1] == 'YTACCATTCAAGCAATGACCTC'
        assert yaml_util.get_arg_as_yaml(wf_inp2["filter_asm_args"])["ctg_len_min"] == 2300

        ## Run w/o existing workflow inputs and with asm_policy spades
        cmd = cmd_base + \
              (" --assembly-policy wgs_spades "
               "--filter-asm-default-ctg-len-min 800")
        check_call(shlex.split(cmd))
        wf_inp3 = "out/gene_extractor.yaml"
        assert os.path.exists(wf_inp3)
        wf_inp3 = yaml_util.load_yaml(wf_inp3)
        print(wf_inp3)
        assert wf_inp3["assembler"] == 'spades'
        assert yaml_util.get_arg_as_yaml(wf_inp3["filter_asm_args"])["ctg_len_min"] == 800

def make_extractor_runner_conf(conda_env,runner_no_retry=True,**kwargs):
    dep_resolver = os.path.abspath("dep_resolver.yaml")
    cmd = "python -m MICGENT.cwl_runner add-dep-resolver --pkg-name ngs-mstb --conda-env {conda_env} --mode w {dep_resolver}".format(
        **locals())
    check_call(shlex.split(cmd))
    kwargs_dfl = dict(
        runner_use_conda=True,
        runner_conda_env="toil",
        dep_resolver=dep_resolver,
        clean_jobstore=True,
        logLevel="DEBUG",
        disableCaching=True
        )

    kwargs_dfl.update(kwargs)
    if runner_no_retry:
        kwargs_dfl.update(dict(retryCount=1,max_wf_tries=1))
    runner_conf = arg_parsing.dump_sig(cwl_runner.run_toil,
                                       out_yaml="cwl_runner.yaml",
                                       kwargs=kwargs_dfl
                                       )
    return runner_conf

def make_extractor_runner_conf_large_run(logLevel="INFO",*l,**kwargs):
    """Use this for large distributed tests to avoid slow donws due to large log transfers"""
    return make_extractor_runner_conf(logLevel=logLevel,*l,**kw)


def path_to_cwl_file(path):
    return {"class": "File","path": path}

@pytest.mark.slow
@pytest.mark.micgent_data
@pytest.mark.large_test_data
@pytest.mark.conda_env_ngs_mstb
@helpers.skip_no_conda_env_toil
def test_clean_reads_small_rsv(request):
    test_data = os.path.abspath(globals()["test_data"])
    data_dir = pjoin(request.config.getoption('--large-test-data'),"rsv/SA2017_sub")
    pkg_data = resources.get_pkg_data_dir("MICGENT")
    cwl_dir = pjoin(pkg_data,"cwl")
    cwl_wf = pjoin(cwl_dir,"clean_reads_qc.cwl")
    conda_env = request.config.getoption('--conda-env-ngs-mstb')
    with helpers.mkchdir("clean_reads_small_rsv"):
        runner_conf = make_extractor_runner_conf(conda_env,runner_no_retry=True) #,batchSystem="Torque")
        inp_reads = [ path_to_cwl_file(os.path.join(data_dir,"SMA1646701_S1_L001_R{}_001.fastq.gz".format(i_read))) \
            for i_read in (1,2) ]
        wf_inp = os.path.abspath("wf_inp.yaml")
        inp = dict(
            SampleID = "SASC_HHC",
            inp_seq1 = inp_reads[0],
            inp_seq2 = inp_reads[1],
            threads = 8
            )
        yaml_util.dump_yaml(inp,wf_inp)
        cmd = ("python -m MICGENT.cwl_runner --config cwl_runner.yaml "
            "run-toil --logLevel DEBUG "
            "{cwl_wf} {wf_inp}").format(**locals())
        check_call(shlex.split(cmd))

@pytest.mark.slow
@pytest.mark.micgent_data
@pytest.mark.large_test_data
@pytest.mark.conda_env_ngs_mstb
@helpers.skip_no_conda_env_toil
def test_ariba_small_sa(request):
    test_data = os.path.abspath(globals()["test_data"])
    conf_dir = pjoin(test_data,"gene_extractor/sa")
    data_dir = pjoin(request.config.getoption('--large-test-data'),"sa/sast2")
    pkg_data = resources.get_pkg_data_dir("MICGENT")
    cwl_dir = pjoin(pkg_data,"cwl")
    cwl_wf = pjoin(cwl_dir,"ariba_run.cwl")
    conda_env = request.config.getoption('--conda-env-ngs-mstb')
    with helpers.mkchdir("ariba_small_sa"):
        runner_conf = make_extractor_runner_conf(conda_env,runner_no_retry=True,batchSystem="singleMachine", leave_workdir=True)
        inp_reads = [ path_to_cwl_file(os.path.join(data_dir,"SMA1828869_S1_L001_R{}_001.fastq.gz".format(i_read))) \
            for i_read in (1,2) ] # SMA1927338_S1_L001_R1_001.fastq.gz
        wf_inp = os.path.abspath("wf_inp.yaml")
        inp = dict(
            reads_1 = inp_reads[0],
            reads_2 = inp_reads[1],
            threads = 8,
            serial = True,
            debug = True,
            SampleID = "SMA1828869",
            prepareref_tgz = path_to_cwl_file(pjoin(conf_dir,"ref_HA.arbref.tar"))
            )
        yaml_util.dump_yaml(inp,wf_inp)
        cmd = ("python -m MICGENT.cwl_runner --config cwl_runner.yaml "
            "run-toil --logLevel DEBUG "
            "{cwl_wf} {wf_inp}").format(**locals())
        check_call(shlex.split(cmd))


@pytest.mark.slow
@pytest.mark.micgent_data
@pytest.mark.large_test_data
@pytest.mark.conda_env_ngs_mstb
@helpers.skip_no_conda_env_toil
def test_ariba_stop_codon_sa(request):
    test_data = os.path.abspath(globals()["test_data"])
    conf_dir = pjoin(test_data,"gene_extractor/sa")
    data_dir = pjoin(request.config.getoption('--large-test-data'),"sa/sast300")
    pkg_data = resources.get_pkg_data_dir("MICGENT")
    cwl_dir = pjoin(pkg_data,"cwl")
    cwl_wf = pjoin(cwl_dir,"ariba_run.cwl")
    conda_env = request.config.getoption('--conda-env-ngs-mstb')
    with helpers.mkchdir("ariba_stop_codon_sa"):
        runner_conf = make_extractor_runner_conf(conda_env,runner_no_retry=True,leave_workdir=True) #batchSystem="Torque",
        inp_reads = [ path_to_cwl_file(os.path.join(data_dir,"SMA807030_S1_L001_R{}_001.fastq.gz".format(i_read))) \
            for i_read in (1,2) ] # SMA1927338_S1_L001_R1_001.fastq.gz
        wf_inp = os.path.abspath("wf_inp.yaml")
        inp = dict(
            reads_1 = inp_reads[0],
            reads_2 = inp_reads[1],
            threads = 8,
            serial = True,
            debug = True,
            SampleID = "SMA807030",
            prepareref_tgz = path_to_cwl_file(pjoin(conf_dir,"ref_HA.arbref.tar"))
            )
        yaml_util.dump_yaml(inp,wf_inp)
        cmd = ("python -m MICGENT.cwl_runner --config cwl_runner.yaml "
            "run-toil --logLevel DEBUG "
            "{cwl_wf} {wf_inp}").format(**locals())
        check_call(shlex.split(cmd))


@pytest.mark.slow
@pytest.mark.micgent_data
@pytest.mark.large_test_data
@pytest.mark.conda_env_ngs_mstb
@helpers.skip_no_conda_env_toil
def test_ariba_no_min_cov_reads_rsv(request):
    test_data = os.path.abspath(globals()["test_data"])
    conf_dir = pjoin(test_data,"gene_extractor/rsv")
    data_dir = pjoin(request.config.getoption('--large-test-data'),"rsv/SA2018")
    pkg_data = resources.get_pkg_data_dir("MICGENT")
    cwl_dir = pjoin(pkg_data,"cwl")
    cwl_wf = pjoin(cwl_dir,"ariba_run.cwl")
    conda_env = request.config.getoption('--conda-env-ngs-mstb')
    with helpers.mkchdir("ariba_run_no_min_cov_reads_rsv"):
        runner_conf = make_extractor_runner_conf(conda_env,runner_no_retry=True,leave_workdir=True)
        inp_reads = [ path_to_cwl_file(os.path.join(data_dir,"SMA1674001_S1_L001_R{}_001.fastq.gz".format(i_read))) \
            for i_read in (1,2) ]
        wf_inp = os.path.abspath("wf_inp.yaml")
        inp = dict(
            reads_1 = inp_reads[0],
            reads_2 = inp_reads[1],
            assembly_cov = 1000000,
            assembly_cov_min = 10,
            assembler = "plugin",
            plugin_asm_options = 'python -m MICGENT.ariba_asm_plugin asm-for-skewed-coverage --extra-args "{deterministic: true}"',
            threads = 8,
            serial = True,
            debug = True,
            SampleID = "SMA1674001",
            prepareref_tgz = path_to_cwl_file(pjoin(conf_dir,"ref_AB_single.trimmed.arbref.tgz"))
            )
        yaml_util.dump_yaml(inp,wf_inp)
        cmd = ("python -m MICGENT.cwl_runner --config cwl_runner.yaml "
            "run-toil --logLevel DEBUG "
            "{cwl_wf} {wf_inp}").format(**locals())
        check_call(shlex.split(cmd))



@pytest.mark.slow
@pytest.mark.micgent_data
@pytest.mark.large_test_data
@pytest.mark.conda_env_ngs_mstb
@helpers.skip_no_conda_env_toil
def test_ariba_one_read_rsv(request):
    test_data = os.path.abspath(globals()["test_data"])
    conf_dir = pjoin(test_data,"gene_extractor/rsv")
    data_dir = pjoin(request.config.getoption('--large-test-data'),"rsv/SA2018_sub")
    pkg_data = resources.get_pkg_data_dir("MICGENT")
    cwl_dir = pjoin(pkg_data,"cwl")
    cwl_wf = pjoin(cwl_dir,"ariba_run.cwl")
    conda_env = request.config.getoption('--conda-env-ngs-mstb')
    with helpers.mkchdir("ariba_run_one_read_rsv"):
        runner_conf = make_extractor_runner_conf(conda_env,runner_no_retry=True,leave_workdir=True)
        inp_reads = [ path_to_cwl_file(os.path.join(data_dir,"SMA1674001_S1_L001_R{}_001.fastq.gz".format(i_read))) \
            for i_read in (1,2) ]
        wf_inp = os.path.abspath("wf_inp.yaml")
        inp = dict(
            reads_1 = inp_reads[0],
            reads_2 = inp_reads[1],
            assembly_cov = 1000000,
            assembly_cov_min = 10,
            assembler = "plugin",
            plugin_asm_options = 'python -m MICGENT.ariba_asm_plugin asm-for-skewed-coverage --extra-args "{deterministic: true}"',
            threads = 8,
            serial = True,
            debug = True,
            SampleID = "SMA1674001",
            prepareref_tgz = path_to_cwl_file(pjoin(conf_dir,"ref_AB_single.trimmed.arbref.tgz"))
            )
        yaml_util.dump_yaml(inp,wf_inp)
        cmd = ("python -m MICGENT.cwl_runner --config cwl_runner.yaml "
            "run-toil --logLevel DEBUG "
            "{cwl_wf} {wf_inp}").format(**locals())
        check_call(shlex.split(cmd))


@pytest.mark.slow
@pytest.mark.micgent_data
@pytest.mark.large_test_data
@pytest.mark.conda_env_ngs_mstb
@helpers.skip_no_conda_env_toil
def test_ariba_zero_reads_rsv(request):
    test_data = os.path.abspath(globals()["test_data"])
    conf_dir = pjoin(test_data,"gene_extractor/rsv")
    data_dir = pjoin(request.config.getoption('--large-test-data'),"rsv/SA2018_sub")
    pkg_data = resources.get_pkg_data_dir("MICGENT")
    cwl_dir = pjoin(pkg_data,"cwl")
    cwl_wf = pjoin(cwl_dir,"ariba_run.cwl")
    conda_env = request.config.getoption('--conda-env-ngs-mstb')
    with helpers.mkchdir("ariba_run_zero_reads_rsv"):
        runner_conf = make_extractor_runner_conf(conda_env,runner_no_retry=True,leave_workdir=True)
        inp_reads = [ path_to_cwl_file(os.path.join(data_dir,"SMA1106665_S1_L001_R{}_001.fastq.gz".format(i_read))) \
            for i_read in (1,2) ]
        wf_inp = os.path.abspath("wf_inp.yaml")
        inp = dict(
            reads_1 = inp_reads[0],
            reads_2 = inp_reads[1],
            assembly_cov = 1000000,
            assembly_cov_min = 10,
            assembler = "plugin",
            plugin_asm_options = 'python -m MICGENT.ariba_asm_plugin asm-for-skewed-coverage --extra-args "{deterministic: true}"',
            threads = 8,
            serial = True,
            debug = True,
            SampleID = "SMA1674001",
            prepareref_tgz = path_to_cwl_file(pjoin(conf_dir,"ref_AB_single.trimmed.arbref.tgz"))
            )
        yaml_util.dump_yaml(inp,wf_inp)
        cmd = ("python -m MICGENT.cwl_runner --config cwl_runner.yaml "
            "run-toil --logLevel DEBUG "
            "{cwl_wf} {wf_inp}").format(**locals())
        check_call(shlex.split(cmd))

@pytest.mark.slow
@pytest.mark.micgent_data
@pytest.mark.large_test_data
@pytest.mark.conda_env_ngs_mstb
@helpers.skip_no_conda_env_toil
def test_ariba_skewed_cov_filter_rsv(request):
    test_data = os.path.abspath(globals()["test_data"])
    conf_dir = pjoin(test_data,"gene_extractor/rsv")
    data_dir = pjoin(request.config.getoption('--large-test-data'),"rsv/SA2018")
    pkg_data = resources.get_pkg_data_dir("MICGENT")
    cwl_dir = pjoin(pkg_data,"cwl")
    cwl_wf = pjoin(cwl_dir,"ariba_run.cwl")
    conda_env = request.config.getoption('--conda-env-ngs-mstb')
    with helpers.mkchdir("ariba_run_skewed_cov_filter_rsv"):
        runner_conf = make_extractor_runner_conf(conda_env,runner_no_retry=True,leave_workdir=True)
        inp_reads = [ path_to_cwl_file(os.path.join(data_dir,"SMA589764_S1_L001_R{}_001.fastq.gz".format(i_read))) \
            for i_read in (1,2) ]
        wf_inp = os.path.abspath("wf_inp.yaml")
        inp = dict(
            reads_1 = inp_reads[0],
            reads_2 = inp_reads[1],
            assembly_cov = 1000000,
            assembly_cov_min = 10,
            assembler = "plugin",
            plugin_asm_options = 'python -m MICGENT.ariba_asm_plugin asm-for-skewed-coverage --extra-args "{deterministic: true}"',
            threads = 8,
            serial = True,
            debug = True,
            SampleID = "SMA589764",
            prepareref_tgz = path_to_cwl_file(pjoin(conf_dir,"ref_AB_single.trimmed.arbref.tgz"))
            )
        yaml_util.dump_yaml(inp,wf_inp)
        cmd = ("python -m MICGENT.cwl_runner --config cwl_runner.yaml "
            "run-toil --logLevel DEBUG "
            "{cwl_wf} {wf_inp}").format(**locals())
        check_call(shlex.split(cmd))


@pytest.mark.slow
@pytest.mark.micgent_data
@pytest.mark.large_test_data
@pytest.mark.conda_env_ngs_mstb
@helpers.skip_no_conda_env_toil
def test_clean_reads_small_sa(request):
    test_data = os.path.abspath(globals()["test_data"])
    data_dir = pjoin(request.config.getoption('--large-test-data'),"sa/sast2")
    pkg_data = resources.get_pkg_data_dir("MICGENT")
    cwl_dir = pjoin(pkg_data,"cwl")
    cwl_wf = pjoin(cwl_dir,"clean_reads_qc.cwl")
    conda_env = request.config.getoption('--conda-env-ngs-mstb')
    with helpers.mkchdir("clean_reads_small_sa"):
        runner_conf = make_extractor_runner_conf(conda_env,runner_no_retry=True,batchSystem="singleMachine")
        inp_reads = [ path_to_cwl_file(os.path.join(data_dir,"SMA1828869_S1_L001_R{}_001.fastq.gz".format(i_read))) \
            for i_read in (1,2) ] # SMA1927338_S1_L001_R1_001.fastq.gz
        wf_inp = os.path.abspath("wf_inp.yaml")
        inp = dict(
            SampleID = "SASMA1828869",
            inp_seq1 = inp_reads[0],
            inp_seq2 = inp_reads[1],
            threads = 8
            )
        yaml_util.dump_yaml(inp,wf_inp)
        cmd = ("python -m MICGENT.cwl_runner --config cwl_runner.yaml "
            "run-toil --logLevel DEBUG "
            "{cwl_wf} {wf_inp}").format(**locals())
        check_call(shlex.split(cmd))


@pytest.mark.slow
@pytest.mark.micgent_data
@pytest.mark.large_test_data
@pytest.mark.conda_env_ngs_mstb
@helpers.skip_no_conda_env_toil
def test_run_extractor_small_rsv(request):
    test_data = os.path.abspath(globals()["test_data"])
    micgent_data = request.config.getoption('--micgent-data')
    conf_dir = pjoin(test_data,"gene_extractor/rsv")
    data_dir = pjoin(request.config.getoption('--large-test-data'),"rsv/SA2017_sub")
    conda_env = request.config.getoption('--conda-env-ngs-mstb')
    with helpers.mkchdir("extractor_small_rsv"):
        runner_conf = make_extractor_runner_conf(conda_env,runner_no_retry=True,batchSystem="singleMachine",clean_jobstore=True,leave_workdir=True)
        cmd = ("python -m MICGENT.gene_extractor run-extraction-wf --micgent-data {micgent_data} "
            "--datadir {data_dir} --prepareref-tgz {conf_dir}/ref_AB_single.trimmed.arbref.tgz "
            "--manifest {conf_dir}/sample_manifest_two.txt --ref-common {conf_dir}/ref_AB_single.fasta "
            "--cwl-runner-name toil --cwl-runner-config cwl_runner.yaml "
            "--cwl-inputs {conf_dir}/cwl_inputs.yaml "
            "--outdir out "
            "--web-dir-out out/web").format(**locals())
        check_call(shlex.split(cmd))
        manifest_out = "out/manifest_out.tsv"
        assert os.path.exists(manifest_out)
        with open(manifest_out,"r") as man:
            lines = man.readlines()
            assert len(lines) == 3, "Expected three lines in output manifest (header and one contig per sample)"
        check_multiqc_cwl_tool_output("out/web/multiqc.html")

@pytest.mark.slow
@pytest.mark.micgent_data
@pytest.mark.large_test_data
@pytest.mark.conda_env_ngs_mstb
@helpers.skip_no_conda_env_toil
def test_run_extractor_no_min_cov_reads_rsv(request):
    test_data = os.path.abspath(globals()["test_data"])
    micgent_data = request.config.getoption('--micgent-data')
    conf_dir = pjoin(test_data,"gene_extractor/rsv")
    data_dir = pjoin(request.config.getoption('--large-test-data'),"rsv/SA2018")
    conda_env = request.config.getoption('--conda-env-ngs-mstb')
    with helpers.mkchdir("extractor_no_min_cov_reads_rsv"):
        runner_conf = make_extractor_runner_conf(conda_env,runner_no_retry=True,leave_workdir=True)
        cmd = ("python -m MICGENT.gene_extractor run-extraction-wf --micgent-data {micgent_data} "
            "--datadir {data_dir} --prepareref-tgz {conf_dir}/ref_AB_single.trimmed.arbref.tgz "
            "--manifest {conf_dir}/sample_manifest_no_min_cov_reads.txt --ref-common {conf_dir}/ref_AB_single.fasta "
            "--cwl-runner-name toil --cwl-runner-config cwl_runner.yaml "
            "--cwl-inputs {conf_dir}/cwl_inputs.yaml "
            "--outdir out "
            "--web-dir-out out/web").format(**locals())
        check_call(shlex.split(cmd))
        manifest_out = "out/manifest_out.tsv"
        assert os.path.exists(manifest_out)
        with open(manifest_out,"r") as man:
            lines = man.readlines()
            assert len(lines) == 1, "Expected only header line in output manifest"
        manifest_out_all = "out/web/manifest_out_all.tsv"
        assert os.path.exists(manifest_out_all)
        man_all = pd.read_table(manifest_out_all)
        assert man_all.shape[0] == 1
        assert 'WARNING:' in str(man_all.Asm_Msg[0])

@pytest.mark.slow
@pytest.mark.micgent_data
@pytest.mark.large_test_data
@pytest.mark.conda_env_ngs_mstb
@helpers.skip_no_conda_env_toil
def test_run_extractor_chim_duplicate_rsv(request):
    """Test that we can drop duplicated matches in chimeric PCR constructs"""
    test_data = os.path.abspath(globals()["test_data"])
    micgent_data = request.config.getoption('--micgent-data')
    conf_dir = pjoin(test_data,"gene_extractor/rsv")
    data_dir = pjoin(request.config.getoption('--large-test-data'),"rsv/SA2018")
    conda_env = request.config.getoption('--conda-env-ngs-mstb')
    with helpers.mkchdir("extractor_chim_duplicate_rsv"):
        runner_conf = make_extractor_runner_conf(conda_env,batchSystem="singleMachine",runner_no_retry=True,leave_workdir=True)
        cmd = ("python -m MICGENT.gene_extractor run-extraction-wf --micgent-data {micgent_data} "
            "--datadir {data_dir} --prepareref-tgz {conf_dir}/ref_AB_single.trimmed.arbref.tgz "
            "--manifest {conf_dir}/sample_manifest_chim_duplicate.txt --ref-common {conf_dir}/ref_AB_single.fasta "
            "--cwl-runner-name toil --cwl-runner-config cwl_runner.yaml "
            "--cwl-inputs {conf_dir}/cwl_inputs.yaml "
            "--outdir out "
            "--web-dir-out out/web").format(**locals())
        check_call(shlex.split(cmd))
        manifest_out = "out/manifest_out.tsv"
        assert os.path.exists(manifest_out)
        man = pd.read_table(manifest_out)
        assert man.shape[0] == 2
        assert man.drop_duplicates(["SampleID"]).shape[0] == 2

@pytest.mark.slow
@pytest.mark.micgent_data
@pytest.mark.large_test_data
@pytest.mark.conda_env_ngs_mstb
@helpers.skip_no_conda_env_toil
def test_run_extractor_skewed_cov_filter_rsv(request):
    """Test that median coverage filter drops sample with very low median base coverage"""
    test_data = os.path.abspath(globals()["test_data"])
    micgent_data = request.config.getoption('--micgent-data')
    conf_dir = pjoin(test_data,"gene_extractor/rsv")
    data_dir = pjoin(request.config.getoption('--large-test-data'),"rsv/SA2018")
    conda_env = request.config.getoption('--conda-env-ngs-mstb')
    with helpers.mkchdir("extractor_skewed_cov_filter_rsv"):
        runner_conf = make_extractor_runner_conf(conda_env,batchSystem="singleMachine",runner_no_retry=True)
        cmd = ("python -m MICGENT.gene_extractor run-extraction-wf --micgent-data {micgent_data} "
            "--datadir {data_dir} --prepareref-tgz {conf_dir}/ref_AB_single.trimmed.arbref.tgz "
            "--manifest {conf_dir}/sample_manifest_skewed_cov_filter.txt --ref-common {conf_dir}/ref_AB_single.fasta "
            "--cwl-runner-name toil --cwl-runner-config cwl_runner.yaml "
            "--cwl-inputs {conf_dir}/cwl_inputs.yaml "
            "--outdir out "
            "--web-dir-out out/web").format(**locals())
        check_call(shlex.split(cmd))
        manifest_out = "out/manifest_out.tsv"
        assert os.path.exists(manifest_out)
        man = pd.read_table(manifest_out)
        assert man.shape[0] == 0

@pytest.mark.slow
@pytest.mark.micgent_data
@pytest.mark.large_test_data
@pytest.mark.conda_env_ngs_mstb
@helpers.skip_no_conda_env_toil
def test_run_extractor_zero_reads_rsv(request):
    test_data = os.path.abspath(globals()["test_data"])
    micgent_data = request.config.getoption('--micgent-data')
    conf_dir = pjoin(test_data,"gene_extractor/rsv")
    data_dir = pjoin(request.config.getoption('--large-test-data'),"rsv/SA2018_sub")
    conda_env = request.config.getoption('--conda-env-ngs-mstb')
    with helpers.mkchdir("extractor_zero_reads_rsv"):
        runner_conf = make_extractor_runner_conf(conda_env,runner_no_retry=True,leave_workdir=True)
        cmd = ("python -m MICGENT.gene_extractor run-extraction-wf --micgent-data {micgent_data} "
            "--datadir {data_dir} --prepareref-tgz {conf_dir}/ref_AB_single.trimmed.arbref.tgz "
            "--manifest {conf_dir}/sample_manifest_zero_reads.txt --ref-common {conf_dir}/ref_AB_single.fasta "
            "--cwl-runner-name toil --cwl-runner-config cwl_runner.yaml "
            "--cwl-inputs {conf_dir}/cwl_inputs.yaml "
            "--outdir out "
            "--web-dir-out out/web").format(**locals())
        check_call(shlex.split(cmd))
        manifest_out = "out/manifest_out.tsv"
        assert os.path.exists(manifest_out)
        with open(manifest_out,"r") as man:
            lines = man.readlines()
            assert len(lines) == 1, "Expected only header line in output manifest"
        manifest_out_all = "out/web/manifest_out_all.tsv"
        assert os.path.exists(manifest_out_all)
        man_all = pd.read_table(manifest_out_all)
        assert man_all.shape[0] == 1
        assert 'WARNING:' in str(man_all.Asm_Msg[0])


def check_multiqc_cwl_tool_output(mqc_html):
    assert not helpers.grep(r"No sample data available to build a report",mqc_html), "Multiqc failed"

@pytest.mark.slow
@pytest.mark.micgent_data
@pytest.mark.large_test_data
@pytest.mark.conda_env_ngs_mstb
@helpers.skip_no_conda_env_toil
def test_run_extractor_small_sa(request):
    test_data = os.path.abspath(globals()["test_data"])
    micgent_data = request.config.getoption('--micgent-data')
    conf_dir = pjoin(test_data,"gene_extractor/sa")
    data_dir = pjoin(request.config.getoption('--large-test-data'),"sa/sast2")
    conda_env = request.config.getoption('--conda-env-ngs-mstb')
    with helpers.mkchdir("extractor_small_sa"):
        runner_conf = make_extractor_runner_conf(conda_env,batchSystem="singleMachine",leave_workdir=True,runner_no_retry=True)
        #runner_conf = make_extractor_runner_conf(conda_env,leave_workdir=True,runner_no_retry=True)
        cmd = ("python -m MICGENT.gene_extractor run-extraction-wf --micgent-data {micgent_data} "
            "--datadir {data_dir} --prepareref-tgz {conf_dir}/ref_HA.arbref.tar "
            "--manifest {conf_dir}/sample_manifest_two.txt --ref-common {conf_dir}/ref_HA.fasta "
            "--cwl-runner-name toil --cwl-runner-config cwl_runner.yaml "
            "--outdir out "
            "--debug "
            "--web-dir-out out/web "
            "--assembly-policy wgs_spades").format(**locals())
        check_call(shlex.split(cmd))
        manifest_out = "out/manifest_out.tsv"
        assert os.path.exists(manifest_out)
        with open(manifest_out,"r") as man:
            lines = man.readlines()
            assert len(lines) == 3, "Expected three lines in output manifest (header and two contigs for one sample)"
        check_multiqc_cwl_tool_output("out/web/multiqc.html")

@pytest.mark.slow
@pytest.mark.micgent_data
@pytest.mark.large_test_data
@pytest.mark.conda_env_ngs_mstb
@helpers.skip_no_conda_env_toil
def test_run_extractor_stop_codon_sa(request):
    test_data = os.path.abspath(globals()["test_data"])
    micgent_data = request.config.getoption('--micgent-data')
    conf_dir = pjoin(test_data,"gene_extractor/sa")
    data_dir = pjoin(request.config.getoption('--large-test-data'),"sa/sast300")
    conda_env = request.config.getoption('--conda-env-ngs-mstb')
    with helpers.mkchdir("extractor_stop_codon_sa"):
        runner_conf = make_extractor_runner_conf(conda_env,batchSystem="singleMachine",leave_workdir=True,runner_no_retry=True)
        #runner_conf = make_extractor_runner_conf(conda_env,leave_workdir=True,runner_no_retry=True)
        cmd = ("python -m MICGENT.gene_extractor run-extraction-wf --micgent-data {micgent_data} "
            "--datadir {data_dir} --prepareref-tgz {conf_dir}/ref_HA.arbref.tar "
            "--manifest {conf_dir}/sample_manifest_sast300_stop_codon.txt --ref-common {conf_dir}/ref_HA.fasta "
            "--cwl-runner-name toil --cwl-runner-config cwl_runner.yaml "
            "--outdir out "
            "--debug "
            "--web-dir-out out/web "
            "--assembly-policy wgs_spades").format(**locals())
        check_call(shlex.split(cmd))
        manifest_out = "out/manifest_out.tsv"
        assert os.path.exists(manifest_out)
        with open(manifest_out,"r") as man:
            lines = man.readlines()
            assert len(lines) == 2, "Expected two lines in output manifest (header and one contig for one sample)"




def run_extractor_from_config(request,
    conf_dir,
    data_dir,
    run_dir,
    man_inp,
    prepareref_tgz,
    ref_common,
    data_root=None,
    man_n_inp=0,
    man_n_repl=10,
    expect_min_rows_manifest_out=0,
    runner_no_retry=False,
    cwl_inputs=None,
    logLevel="INFO",
    leave_workdir=None):
    """Generic test function for config-driven deterministic test.
    :param expect_min_rows_manifest_out: Check that at least that many rows are returned per non-replicated output manifest.
    :param run_dir: relative to the callers cwd. Output will be in run_dir/out.

    All other input paths should be relative to respective top data directories of gene_extractor tests.
    """
    test_data = os.path.abspath(globals()["test_data"])
    micgent_data = request.config.getoption('--micgent-data')
    extra_config = request.config.getoption('--extra-config')
    if extra_config:
        extra_config = yaml_util.load_yaml(extra_config)
    else:
        extra_config = {}
    conf_dir = pjoin(test_data,"gene_extractor",conf_dir)
    if not cwl_inputs:
        cwl_inputs = "cwl_inputs.yaml"
    cwl_inputs = pjoin(conf_dir,cwl_inputs)
    if data_root is None:
        data_root = request.config.getoption('--large-test-data')
    data_dir = pjoin(data_root,data_dir)
    conda_env = request.config.getoption('--conda-env-ngs-mstb')
    with helpers.mkchdir(run_dir):
        man_inp = "{conf_dir}/{man_inp}".format(**locals())
        man_rep_df = ge_bench.replicate_sample_man(man_inp,n_inp=man_n_inp,n_repl=man_n_repl)
        man_rep = os.path.join(os.getcwd(),"sample_manifest_rep.txt")
        man_rep_df.to_csv(man_rep, index=False, sep="\t")
        runner_kw = extra_config.get("cwl_runner",{}).get("run_toil",{})
        if logLevel:
            runner_kw["logLevel"] = logLevel
        if leave_workdir is not None:
            runner_kw["leave_workdir"] = leave_workdir
        runner_conf = make_extractor_runner_conf(conda_env,
            runner_no_retry=runner_no_retry,
            **runner_kw)
        cmd = ("python -m MICGENT.gene_extractor run-extraction-wf --micgent-data {micgent_data} "
            "--datadir {data_dir} --prepareref-tgz {conf_dir}/{prepareref_tgz} "
            "--manifest {man_rep} --ref-common {conf_dir}/{ref_common} "
            "--cwl-runner-name toil --cwl-runner-config cwl_runner.yaml "
            "--cwl-inputs {cwl_inputs} "
            "--deterministic "
            "--outdir out "
            "--web-dir-out out/web").format(**locals())
        check_call(shlex.split(cmd))
        if man_rep_df.shape[0] > 0:
            check_multiqc_cwl_tool_output("out/web/multiqc.html")
        manifest_out = "out/manifest_out.tsv"
        assert os.path.exists(manifest_out)
        man_out = pd.read_table(manifest_out)
        assert man_out.shape[0] >= expect_min_rows_manifest_out * man_n_repl

def run_extractor_determinism_multi_run(run_dir,n_runs=2,*l,**kw):
    """Determinism test by replicating within run and doing multiple runsself.
    """
    run_dir_test=run_dir
    #with helpers.mkchdir(run_dir_test):
    run_out_dirs = []
    for i_run in range(1,n_runs+1):
        run_dir = pjoin(run_dir_test,"r{}".format(i_run))
        run_extractor_from_config(run_dir=run_dir,*l,**kw)
        out_dir = pjoin(run_dir,"out")
        ge_bench.check_within_run_replicate_equality(out_dir,no_assert=False)
        run_out_dirs.append(out_dir)

    ge_bench.check_runs_for_equality(run_out_dirs,no_assert=False,use_existing_sigs=True)


@pytest.mark.slow
@pytest.mark.micgent_data
@pytest.mark.large_test_data
@pytest.mark.conda_env_ngs_mstb
@helpers.skip_no_conda_env_toil
def test_run_extractor_determinism_rsv_1(request):
    run_dir="extractor_determinism_rsv_1"
    run_extractor_from_config(request=request,
        conf_dir="rsv",
        data_dir="rsv/SA2017_sub",
        run_dir=run_dir,
        man_inp="sample_manifest_two.txt",
        prepareref_tgz="ref_AB_single.trimmed.arbref.tgz",
        ref_common="ref_AB_single.fasta",
        man_n_inp=0,
        man_n_repl=10,
        leave_workdir=True)
    ge_bench.check_within_run_replicate_equality(pjoin(run_dir,"out"),no_assert=False)

@pytest.mark.slow
@pytest.mark.micgent_data
@pytest.mark.large_test_data
@pytest.mark.conda_env_ngs_mstb
@helpers.skip_no_conda_env_toil
def test_run_extractor_determinism_rsv_complete_1(request):
    run_dir="extractor_determinism_rsv_complete_1"
    run_extractor_from_config(request=request,
        conf_dir="rsv",
        cwl_inputs="cwl_inputs_complete.yaml",
        data_dir="rsv/INF2019",
        run_dir=run_dir,
        man_inp="sample_manifest_complete_one.txt",
        prepareref_tgz="ref_AB_singleWG_RSV.arbref.tgz",
        ref_common="ref_AB_singleWG_RSV.fasta",
        man_n_inp=0,
        man_n_repl=2,
        expect_min_rows_manifest_out=1)
    ge_bench.check_within_run_replicate_equality(pjoin(run_dir,"out"),no_assert=False)

@pytest.mark.slow
@pytest.mark.micgent_data
@pytest.mark.large_test_data
@pytest.mark.conda_env_ngs_mstb
@helpers.skip_no_conda_env_toil
def test_run_extractor_determinism_sa_1(request):
    run_dir="extractor_determinism_sa_1"
    run_extractor_from_config(request=request,
        conf_dir="sa",
        data_dir="sa/sast2",
        run_dir=run_dir,
        man_inp="sample_manifest_two.txt",
        prepareref_tgz="ref_HA.arbref.tar",
        ref_common="ref_HA.fasta",
        man_n_inp=0,
        man_n_repl=10)
    ge_bench.check_within_run_replicate_equality(pjoin(run_dir,"out"),no_assert=False)

@pytest.mark.slow
@pytest.mark.micgent_data
@pytest.mark.large_test_data
@pytest.mark.conda_env_ngs_mstb
@helpers.skip_no_conda_env_toil
def test_run_extractor_determinism_sa_2(request):
    run_dir="extractor_determinism_sa_2"
    run_extractor_from_config(request=request,
        conf_dir="sa",
        data_dir="sa/sast300",
        run_dir=run_dir,
        man_inp="sample_manifest_sast300_tough.txt",
        prepareref_tgz="ref_HA.arbref.tar",
        ref_common="ref_HA.fasta",
        man_n_inp=0,
        man_n_repl=5,
        leave_workdir=True)
    ge_bench.check_within_run_replicate_equality(pjoin(run_dir,"out"),no_assert=False)

@pytest.mark.slow
@pytest.mark.micgent_data
@pytest.mark.large_test_data
@pytest.mark.conda_env_ngs_mstb
@helpers.skip_no_conda_env_toil
def test_run_extractor_determinism_sa_fermi_2(request):
    run_dir="extractor_determinism_sa_fermi_2"
    run_extractor_from_config(request=request,
        conf_dir="sa",
        cwl_inputs="cwl_inputs_fermilite.yaml",
        data_dir="sa/sast300",
        run_dir=run_dir,
        man_inp="sample_manifest_sast300_tough.txt",
        prepareref_tgz="ref_HA.arbref.tar",
        ref_common="ref_HA.fasta",
        man_n_inp=0,
        man_n_repl=5,
        leave_workdir=True)
    ge_bench.check_within_run_replicate_equality(pjoin(run_dir,"out"),no_assert=False)

@pytest.mark.slow
@pytest.mark.micgent_data
@pytest.mark.large_test_data
@pytest.mark.conda_env_ngs_mstb
@helpers.skip_no_conda_env_toil
def test_run_extractor_determinism_multi_run_sa_1(request):
    run_dir="extractor_determinism_multi_run_sa_1"
    run_extractor_determinism_multi_run(
        n_runs=2,
        request=request,
        conf_dir="sa",
        data_dir="sa/sast300",
        run_dir=run_dir,
        man_inp="sample_manifest_sast300_tough.txt",
        prepareref_tgz="ref_HA.arbref.tar",
        ref_common="ref_HA.fasta",
        man_n_inp=0,
        man_n_repl=2)

@pytest.mark.slow
@pytest.mark.micgent_data
@pytest.mark.large_test_data
@pytest.mark.conda_env_ngs_mstb
@helpers.skip_no_conda_env_toil
def test_run_extractor_determinism_multi_run_sa_2(request):
    run_dir="extractor_determinism_multi_run_sa_2"
    run_extractor_determinism_multi_run(
        n_runs=2,
        request=request,
        conf_dir="sa",
        cwl_inputs="cwl_inputs_fermilite.yaml",
        data_dir="sa/sast300",
        run_dir=run_dir,
        man_inp="sample_manifest_sast300_tough.txt",
        prepareref_tgz="ref_HA.arbref.tar",
        ref_common="ref_HA.fasta",
        man_n_inp=0,
        man_n_repl=2)

@pytest.mark.slow
@pytest.mark.micgent_data
@pytest.mark.large_test_data
@pytest.mark.conda_env_ngs_mstb
@helpers.skip_no_conda_env_toil
def test_run_extractor_determinism_multi_run_rsv_1(request):
    run_dir="extractor_determinism_multi_run_rsv_1"
    run_extractor_determinism_multi_run(
        n_runs=2,
        request=request,
        conf_dir="rsv",
        data_dir="rsv/SA2017_sub",
        run_dir=run_dir,
        man_inp="sample_manifest_two.txt",
        prepareref_tgz="ref_AB_single.trimmed.arbref.tgz",
        ref_common="ref_AB_single.fasta",
        man_n_inp=0,
        man_n_repl=2)

@pytest.mark.slow
@pytest.mark.micgent_data
@pytest.mark.huge_test_data
@pytest.mark.conda_env_ngs_mstb
@helpers.skip_no_conda_env_toil
def test_run_extractor_determinism_multi_run_huge_rsv_2(request):
    data_root = request.config.getoption('--huge-test-data')
    run_dir="extractor_determinism_multi_run_huge_rsv_2"
    run_extractor_determinism_multi_run(
        data_root=data_root,
        n_runs=2,
        request=request,
        conf_dir="rsv",
        data_dir="rsv/SA2017",
        run_dir=run_dir,
        man_inp="sample_manifest_repl_hard.txt",
        prepareref_tgz="ref_AB_single.trimmed.arbref.tgz",
        ref_common="ref_AB_single.fasta",
        man_n_inp=0,
        man_n_repl=100)

@pytest.mark.slow
@pytest.mark.micgent_data
@pytest.mark.huge_test_data
@pytest.mark.conda_env_ngs_mstb
@helpers.skip_no_conda_env_toil
def test_run_extractor_determinism_multi_run_huge_rsv_1(request):
    data_root = request.config.getoption('--huge-test-data')
    run_dir="extractor_determinism_multi_run_huge_rsv_1"
    run_extractor_determinism_multi_run(
        data_root=data_root,
        n_runs=2,
        request=request,
        conf_dir="rsv",
        data_dir="rsv/SA2017",
        run_dir=run_dir,
        man_inp="sample_manifest.txt",
        prepareref_tgz="ref_AB_single.trimmed.arbref.tgz",
        ref_common="ref_AB_single.fasta",
        man_n_inp=0,
        man_n_repl=2,
        expect_min_rows_manifest_out=145)

@pytest.mark.slow
@pytest.mark.micgent_data
@pytest.mark.huge_test_data
@pytest.mark.conda_env_ngs_mstb
@helpers.skip_no_conda_env_toil
def test_run_extractor_determinism_multi_run_huge_sa_1(request):
    data_root = request.config.getoption('--huge-test-data')
    run_dir="extractor_determinism_multi_run_huge_sa_1"
    run_extractor_determinism_multi_run(
        data_root=data_root,
        n_runs=2,
        request=request,
        conf_dir="sa",
        data_dir="sa/sast300",
        run_dir=run_dir,
        man_inp="sample_manifest.txt",
        prepareref_tgz="ref_HA.arbref.tar",
        ref_common="ref_HA.fasta",
        man_n_inp=0,
        man_n_repl=2)

@pytest.mark.slow
@pytest.mark.micgent_data
@pytest.mark.huge_test_data
@pytest.mark.conda_env_ngs_mstb
@helpers.skip_no_conda_env_toil
def test_run_extractor_determinism_multi_run_huge_sa_2(request):
    data_root = request.config.getoption('--huge-test-data')
    run_dir="extractor_determinism_multi_run_huge_sa_2"
    run_extractor_determinism_multi_run(
        data_root=data_root,
        n_runs=2,
        request=request,
        conf_dir="sa",
        cwl_inputs="cwl_inputs_fermilite.yaml",
        data_dir="sa/sast300",
        run_dir=run_dir,
        man_inp="sample_manifest.txt",
        prepareref_tgz="ref_HA.arbref.tar",
        ref_common="ref_HA.fasta",
        man_n_inp=0,
        man_n_repl=2,
        expect_min_rows_manifest_out=290)
