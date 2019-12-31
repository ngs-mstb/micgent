from __future__ import absolute_import
from __future__ import unicode_literals
from __future__ import print_function
from __future__ import division
from future import standard_library
standard_library.install_aliases()
from builtins import str
from builtins import *
from . import util
from six.moves import shlex_quote

import subprocess
from subprocess import check_call
import os
import shlex
import shutil
import sys

from MGT.Logging import *

def cmd_to_str(cmd):
    if util.is_sequence(cmd):
        cmd = " ".join(shlex_quote(str(_)) for _ in cmd)
    return cmd

def get_conda_root(env=None):
    conda_root = None
    conda_exe = shutil.which("conda")
    if conda_exe:
        conda_root = os.path.dirname(os.path.dirname(conda_exe))
        log.debug("Obtained conda_root from conda_exe: {}".format(conda_exe))
    else:
        if env is None:
            env = os.environ
        conda_python = env.get("CONDA_PYTHON_EXE",None)
        if conda_python:
            conda_root = os.path.dirname(os.path.dirname(conda_python))
            log.debug("Obtained conda_root from conda_python: {}".format(conda_python))
    log.debug("Obtained conda_root as: {}".format(conda_root))
    return conda_root

def conda_activate_cmd(conda_root=None,conda_env=None,conda_activate="conda activate",env=None,require_exist=False):
    """
    Generate BASH command line to activate Conda environment from arbitrary Conda installation.

    :param use_conda: Execute the command in a Conda environment. This triggers the processing of
    other Conda-related arguments, including a search for a current Conda environment if necessary.
    :param conda_root: Installation directory of Conda, must exist if provided. If use_conda is True,
    and conda_root is not provided, Conda should be already initialized in the calling process so that
    this function could detect the active Conda installation and activate it for the subprocess. The latter is
    necessary if conda_activate=="conda activate" (new Conda way of activating the environments), because somehow
    the state of "conda activate" does not propagate fully into a sub-shell.
    :param conda_env: Conda environment to activate before the command is executed
    :param conda_activate: ["conda activate"|"source activate"]. If conda_env is not provided, but use_conda is True,
    the root environment (called "base" in the latest versions of Conda) will be activated with this command.
    :param env: dict to extract environment variables from; None means os.environ
    :param require_exist: If True, assert that Conda profile file exists locally if defined
    :return: One-line command string
    """
    cmd_pref = ""
    if not conda_root:
        ## this would be defined after sourcing conda_root/etc/profile.d/conda.sh or activating the environment
        conda_root = get_conda_root(env=env)
    if conda_root:
        if require_exist:
            assert os.path.isdir(conda_root), "conda_root directory does not exist or is not a directory: {}".format(
            conda_root)
        conda_profile = os.path.join(conda_root, "etc/profile.d/conda.sh")
        if os.path.exists(conda_profile):
            cmd_pref = ". {}; ".format(conda_profile)
        ## we need to manipulate PATH in the command string and not in env_arg because
        ## the shell if bash for conda use, and can source the profiles resetting the PATH
        else:
            cmd_pref = 'export PATH="{}:$PATH"; '.format(os.path.join(conda_root, "bin"))
    ## both 'conda activate' and 'source activate' work fine with no argument, activating
    ## 'base' or 'root' environment respectively
    if not conda_env:
        conda_env = ""
    cmd_pref += "{} {}; ".format(conda_activate, conda_env)
    log.debug("Conda activation command: {}".format(cmd_pref))
    return cmd_pref


step_default_shell = "/bin/bash --noprofile --norc"

def run_step(cmd,step_dir=".",descr=None,stdout=None,stderr=None,args=None,
             shell=False,env=None,env_merge=True,rm_step_dir=False,rm_step_dir_before=False,
             shell_exe=step_default_shell,
             shell_strict=True,
             use_conda=False,
             conda_root=None,
             conda_env=None,
             conda_activate="conda activate",
             dry_run=False):
    """Wrapper method around subprocess.check_call to run a command in a subprocess.
    Can optionally activate Conda environment from an arbitrary Conda installation directory
    as well as add / or replace environment variables. If executing through a shell, uses
    /bin/sh by default in order to minimize side-effects from sourcing user profiles by Bash.
    This function constructs the arguments for subprocess.check_call, and then calls that method.

    :param cmd: Command to run, either a list or a shell-quoted string
    :param step_dir: Will try to create this directory (with intermediates) if does not exist already, and
    change into it when running the command. No attempt is made to adjust relative paths in the command if
    such are present
    :param rm_step_dir: Will try to delete the step_dir upon command completion unless step_dir is the
    current working directory when this method is called. In combination with step_dir, this allows executing
    the command in a temporary directory that is wiped-out afterwards.
    :param use_conda: Execute the command in a Conda environment. This triggers the processing of
    other Conda-related arguments, including a search for a current Conda environment if necessary.
    :param conda_root: Installation directory of Conda, must exist if provided. If use_conda is True,
    and conda_root is not provided, Conda should be already initialized in the calling process so that
    this function could detect the active Conda installation and activate it for the subprocess. The latter is
    necessary if conda_activate=="conda activate" (new Conda way of activating the environments), because somehow
    the state of "conda activate" does not propagate fully into a sub-shell.
    :param conda_env: Conda environment to activate before the command is executed
    :param conda_activate: ["conda activate"|"source activate"]. If conda_env is not provided, but use_conda is True,
    the root environment (called "base" in the latest versions of Conda) will be activated with this command.
    :param dry_run: Do not execute the command, only return the constructed check_call arguments

    :return: dict with the constructed check_call arguments
    """
    ##TODO: consider replacing this method with a call to cwltool with the constructed on-the-fly command tool and
    ##dependency resultion files. This will allow using all dependency resolution mechanisms supported by cwltool,
    ##and not just Conda
    stdout_close = False
    stderr_close = False
    if stdout is not None and util.is_string(stdout):
        if stdout != "stdout":
            stdout = open(stdout,"w")
            stdout_close = True
        else:
            stdout = None
    if stderr is not None and util.is_string(stderr):
        if stderr != "stderr":
            stderr = open(stderr,"w")
            stderr_close = True
        else:
            stderr = None
    try:
        if descr is None:
            descr = str(cmd)
        if args is None:
            args = {}
        if use_conda:
            shell = True
            shell_exe = step_default_shell
        if stdout:
            print("## {}".format(descr), file=stdout)
            stdout.flush()
        if env is not None:
            if env_merge:
                env_arg = os.environ.copy()
                env_arg.update(env)
            else:
                env_arg = env.copy()
        else:
            env_arg = None
        if shell or args:
            cmd = cmd_to_str(cmd)
            cmd = cmd.format(**args)
        if not shell:
            ## presence of args might have already converted cmd to str
            if util.is_string(cmd):
                cmd = shlex.split(cmd)
        else:
            ## build command prefix
            cmd_pref = ""

            if use_conda:
                cmd_pref += conda_activate_cmd(conda_root=conda_root,conda_env=conda_env,conda_activate=conda_activate,require_exist=True)

            ## cmd is already a str
            if shell_strict:
                cmd = "set -euo pipefail; " + cmd
            cmd = cmd_pref + cmd
            ## By default, we want /bin/sh, not /bin/bash, in order
            ## to avoid sourcing bashrc and potentially screwing the environment.
            if not shell_exe:
                shell_exe = step_default_shell
            if util.is_string(shell_exe):
                shell_exe = shlex.split(shell_exe)
            cmd = shell_exe +["-c",cmd]
        if not dry_run:
            if not os.path.exists(step_dir):
                os.makedirs(step_dir,exist_ok=True)
            else:
                if rm_step_dir_before and not os.path.samefile(os.getcwd(), step_dir):
                    shutil.rmtree(step_dir, ignore_errors=True)
                    os.makedirs(step_dir, exist_ok=True)
            check_call(cmd, stdout=stdout, stderr=stderr,cwd=step_dir, env=env_arg)
            if rm_step_dir and not os.path.samefile(os.getcwd(),step_dir):
                shutil.rmtree(step_dir,ignore_errors=True)
        log.info("run_step command: {}".format(cmd))
    finally:
        if stdout_close:
            stdout.close()
        if stderr_close:
            stderr.close()
        #log.debug(subprocess.run("vmstat",stdout=subprocess.PIPE,encoding="utf8").stdout)
            
    return dict(cmd=cmd,env=env_arg,cwd=step_dir)


bbtools_bool = lambda x: "true" if  x else "false"

bbtools_kw_translate = lambda x: "in" if  x=="inp" else "in2" if x=="inp2" else x

def bbtools_abs_file(fname,specials=(r"stdin\.",r"stdout\.")):
    if fname:
        if util.re_match_any(specials,fname):
            return fname
        else:
            return os.path.abspath(fname)
    return fname

def bbtools_abs_adapter_file(fname,specials=('adapters','artifacts','phix','lambda','pjet','mtst','kapa')):
    return bbtools_abs_file(fname,specials=specials)

def bbtools_cmd(exe=None,
    threads=1,
    ram=4*1024,
    rams=2*1024,
    overwrite=True,
    eoom=True,
    **kw):
    kw = kw.copy()
    s = []
    if exe:
        s += [exe]
    if ram:
        s += ["-Xmx{}m".format(ram)]
    if rams:
        if ram and (rams>ram):
            rams = ram
        s += ["-Xms{}m".format(rams)]
    if eoom:
        s += ["-eoom"]
    kw["threads"] = threads
    kw["overwrite"] = overwrite
    for k,v in kw.items():
        if v is not None:
            if util.is_bool(v):
                v = bbtools_bool(v)            
            s += ["{}={}".format(bbtools_kw_translate(k),v)]
    return s

def bbduk_cmd(exe="bbduk.sh",statscolumns=5,ordered=True,tbo="t",copyundefined="t",k=23,mink=2,
    hdist=1,hdist2=0,**kw):
    args = locals().copy()
    del args["kw"]
    args.update(kw)
    return bbtools_cmd(**args)

def reformat_subsample_stdinp(exe="reformat.sh",sampleseed=123,
    inp="stdin.fq",
    threads=1, ram=1024, interleaved=True,
    samplerate=0.1,
    out=None,out2=None,
    **kw):
    """Using Unix `tee`, divert a subset of reads from stdin into a file for QC.
    """
    args = locals().copy()
    del args["kw"]
    args.update(kw)    
    assert inp.startswith("stdin."), "Need standard input with the BBTools format specifier"
    ## If needed, tee can redirect to multiple subcommands: tee >(cmd1) >(cmd2) | ...
    cmd = "tee >(" + cmd_to_str(bbtools_cmd(**args)) + ")"
    return cmd


def clean_reads(inp_reads="stdin.fq.gz",
                out_reads="stdout.fq.gz",
                inp_reads2=None,
                out_reads2=None,
                base="cln",
                workdir="cln",
                adapter_file=None,
                primer_literals=None,
                spikes_file=None,
                minlen=13,
                maq=None,
                trimq=None,
                qtrim=None,
                deterministic=False,
                clumpify=False,
                filter_spikes=False,
                threads=1,
                gc_threads=None,
                ram=4*1024,
                inp_interleaved=False,
                out_stats=None,
                out_qc_before_reads=None,
                out_qc_before_reads2=None,
                out_qc_after_reads=None,
                out_qc_after_reads2=None,                
                out_qc_samplerate=0.1,
                dry_run=False,
                stdout=None,
                stderr=None
                ):
    """Clean reads - trim adapters, primers, remove spikes, optical duplicates.
    All operations are optional.
    Returned paths are absolute"""
    ## WARNING: if later implementing here subsampling options, make sure to
    ## sort **before** the subsampling if deterministic==True
    if not stdout:
        stdout = os.path.abspath("{}.step.stdout".format(base))
    if not stderr:
        stderr = os.path.abspath("{}.step.stderr".format(base))
    if not gc_threads:
        gc_threads = threads
    ## Convert input paths to absolute paths
    inp_reads = bbtools_abs_file(inp_reads)
    inp_reads2 = bbtools_abs_file(inp_reads2)    
    if inp_reads2:
        assert inp_interleaved == False, "Cannot have two input files in interleaved format"
    inp_ilv = inp_interleaved
    adapter_file = bbtools_abs_adapter_file(adapter_file)
    spikes_file = bbtools_abs_file(spikes_file)
    out_reads = bbtools_abs_file(out_reads)
    out_reads2 = bbtools_abs_file(out_reads2)
    
    out_qc_before_reads = bbtools_abs_file(out_qc_before_reads)
    out_qc_before_reads2 = bbtools_abs_file(out_qc_before_reads2)
    out_qc_after_reads = bbtools_abs_file(out_qc_after_reads)
    out_qc_after_reads2 = bbtools_abs_file(out_qc_after_reads2)

    out_stats = bbtools_abs_file(out_stats)

    with util.chdir(workdir,create=True,to_abs=True) as dh_work:
        _read_names = dh_work.make_paired_files_func(base="read",ext="fastq.gz",n=2)
        _read_name = dh_work.make_file_func(base="read",ext="fastq")
        _read_name_gz = dh_work.make_file_func(base="read",ext="fastq.gz")
        _file_name = dh_work.make_file
        stats_files = []
        cmd = cmd_to_str(bbtools_cmd(exe="reformat.sh", sampleseed=123, inp=inp_reads, inp2=inp_reads2,
            out="stdout.fq", interleaved=inp_ilv, ram=ram))
        if inp_reads2:
            inp_ilv = True

        if out_qc_before_reads and out_qc_samplerate > 0:
            cmd += " | " + reformat_subsample_stdinp(interleaved=inp_ilv,
                samplerate=out_qc_samplerate,
                out=out_qc_before_reads,out2=out_qc_before_reads2)
        if clumpify:
            out_clump = _read_name_gz("clump")
            ## TODO: expose those dedup parameters that have to be changed by the user according to the instrument model
            ## https://www.biostars.org/p/225338/
            ## For dedup, as the author advises, it is best to keep subs param at a non-zero value, otherwise will enrich dataset
            ## with reads containing errors.            
            cmd += " | (" + cmd_to_str(bbtools_cmd(exe="clumpify.sh", seed=123, 
                dedupe=True, optical=True, 
                threads=threads, inp="stdin.fq", interleaved=inp_ilv, ram=ram, out=out_clump))
            cmd += " && zcat {out_clump}) ".format(**locals())
        
        if filter_spikes:
            if not spikes_file:
                spikes_file = "phix"
            spikes_stats = _file_name("bbduk_clump","stats")
            stats_files.append(spikes_stats)
            cmd += " | " + cmd_to_str(bbduk_cmd(ref=spikes_file,ktrim=False,k=31,mink=-1,hdist=1,tbo=False,
                stats=spikes_stats,
                threads=threads,interleaved=inp_ilv, ram=ram, inp="stdin.fastq",out="stdout.fq"))
        if not adapter_file:
            adapter_file = "adapters"
        ada_r_stats = _file_name("bbduk_ada_r","stats")
        stats_files.append(ada_r_stats)
        
        do_qc_sub_after = (out_qc_after_reads and out_qc_samplerate > 0)
        # default k, mink and hdist are taken from bbduk author's post on adapter trimming
        # https://www.biostars.org/p/155165/#268947
        # Except that we set mink=2 in order to catch partial read-through adapters
        # at the expense of losing some small amount of legitimate end fragments (which
        # are usually low quality anyway). hdist2 is set to 0.
        
        bbd_kw_last = dict(minlen=minlen)
        if maq is not None:
            bbd_kw_last["maq"] = maq
        if trimq is not None:
            bbd_kw_last["trimq"] = trimq
        if qtrim is not None:
            bbd_kw_last["qtrim"] = qtrim
        out_kw_last = dict(out=out_reads)
        if out_reads2:
            out_kw_last["out2"] = out_reads2
        if not (deterministic or do_qc_sub_after):
            bbd_kw_last.update(out_kw_last)
        else:
            bbd_kw_last["out"] = "stdout.fq"

        if not primer_literals:
            bb_kw = bbd_kw_last
        else:
            bb_kw = dict(out="stdout.fq")
        cmd += " | " + cmd_to_str(bbduk_cmd(ref=adapter_file,ktrim="r",rcomp="t",hdist=1,tpe="t",ftm=5,
            stats=ada_r_stats,
            threads=threads,interleaved=inp_ilv, ram=ram, inp="stdin.fastq",
            **bb_kw))
        if primer_literals:
            if util.is_sequence(primer_literals):
                primer_literals = ",".join(primer_literals)
            pr_l_stats = _file_name("bbduk_pr_l","stats")
            stats_files.append(pr_l_stats)
            pr_lens = [ len(_) for _ in primer_literals.split(",") ]
            pr_k = min(pr_lens)
            pr_restrictleft = 2*max(pr_lens)
            cmd += " | " + cmd_to_str(bbduk_cmd(literal=primer_literals,ktrim="l",rcomp="t",mm="f",k=pr_k,mink=2,
                restrictleft=pr_restrictleft,
                stats=pr_l_stats,
                threads=threads,interleaved=inp_ilv, ram=ram, inp="stdin.fastq",out="stdout.fq"))
            pr_r_stats = _file_name("bbduk_pr_r","stats")
            stats_files.append(pr_r_stats)            
            cmd += " | " + cmd_to_str(bbduk_cmd(literal=primer_literals,ktrim="r",rcomp="t",mm="f",k=pr_k,mink=2,
                stats=pr_r_stats,
                threads=threads,interleaved=inp_ilv, ram=ram, inp="stdin.fastq",
                **bbd_kw_last))
        
        if do_qc_sub_after:
            cmd += " | " + reformat_subsample_stdinp(interleaved=inp_ilv,
                samplerate=out_qc_samplerate,
                out=out_qc_after_reads,out2=out_qc_after_reads2)
            
            if not deterministic:
                cmd += " | " + cmd_to_str(bbtools_cmd(exe="reformat.sh", sampleseed=123, inp="stdin.fastq",
                    interleaved=inp_ilv, ram=1024,
                    **out_kw_last))

        if deterministic:
            ## sortbyname.sh cannot split interleaved reads, so we need to add a reformat.sh step in such case
            ## sortbyname.sh <= v.38.* hangs with an exception when using temporary sorting files
            ## on interleaved inputs. We need to give it enough RAM and disable temp files.
            ## Temp files seem to work OK on a local FS on Mac, but fail on NFS on Linux
            sort_ram = max(1024*4,ram)
            if not inp_ilv:
                cmd += " | " + cmd_to_str(bbtools_cmd(exe="sortbyname.sh",inp="stdin.fastq",
                    threads=threads,interleaved=inp_ilv, allowtemp=False, ram=sort_ram,
                    **out_kw_last))
            else:
                cmd += " | " + cmd_to_str(bbtools_cmd(exe="sortbyname.sh",inp="stdin.fq",
                    threads=threads,interleaved=inp_ilv, allowtemp=False, ram=sort_ram,
                    out="stdout.fq"))
                cmd += " | " + cmd_to_str(bbtools_cmd(exe="reformat.sh", sampleseed=123, inp="stdin.fq",
                    interleaved=inp_ilv, ram=1024,
                    **out_kw_last))
        run_step(cmd,
                 descr="Clean reads",
                 stdout=stdout, stderr=stderr, args=locals(),
                 shell=True,
                 rm_step_dir=False,
                 ## it is citical that we wipe out any RAM settings in _JAVA_OPTIONS,
                 ## so that the ram arguments in the individual tool CLIs would take effect
                 env=dict(_JAVA_OPTIONS="-XX:ParallelGCThreads={gc_threads} -XX:+UseParallelOldGC".format(gc_threads=gc_threads)),
                 dry_run=dry_run)
        if not dry_run and out_stats:
            util.cat_files(stats_files,out_stats)

def replace_fasta_ids(fasta_inp,fasta_out,prefix,key1="1",key2="1"):
    key1 = str(key1)
    key2 = str(key2)
    assert not (("." in key1) or ("." in key2)), "Componenst of contig ID should not contain dots"
    header_tpl = ">{}.l{}.c{}.".format(prefix,key1,key2)+"ctg.{}\n"
    with open(fasta_inp, "r") as inp, open(fasta_out, "w") as out:
        n_cont = 0
        for line in inp:
            if line.startswith(">"):
                n_cont += 1
                line = header_tpl.format(n_cont)
            out.write(line)
    return n_cont

def sort_reads(reads,out_reads,stdout=None,stderr=None,threads=1,base="rd_srt"):
    for i_rd, inp_rd, out_rd in zip(range(1, len(reads) + 1),
                                    reads, out_reads):
        run_step("seqkit sort --threads {threads} -o '{out_rd}' '{inp_rd}'",
                 step_dir="{}_{}".format(base,i_rd),
                 descr="Sort reads, read file No: {}".format(i_rd),
                 stdout=stdout, stderr=stderr, args=locals())


def bwa_mem_long(reads,ref,stdout=None,stderr=None,threads=1,base="bwa_mem",workdir="."):
    """Map long sequences with bwa mem.
    Returned paths are absolute"""
    ref = os.path.abspath(ref)
    reads = os.path.abspath(reads)
    with util.chdir(workdir,create=True,to_abs=True) as dh_work:
        _file_name = dh_work.make_file
        bam = _file_name(base, "bam")
        step_dir = "{}.tmp".format(base)
        ## 0x800 are supplementary alignments
        cmd =  """\
        ln -s {ref} ref.fa
        bwa index ref.fa && \\
        (bwa mem -t {threads} \\
         -w 1000 ref.fa {reads} \\
        | samtools view -@ {threads} -bSu - \\
        | samtools view -b -F 0x800 - \\
        | samtools sort -@ {threads} -T {reads}.bwa.tmp -o {bam} -) && \\
        samtools index {bam} && \\
        rm -rf {reads}.bwa.tmp
        """
        run_step(cmd,
                 step_dir=step_dir,
                 descr="Map long sequences with BWA MEM",
                 stdout=stdout, stderr=stderr, args=locals(),
                 shell=True,
                 rm_step_dir=True)
    return dict(bam=bam)

def minimap2_contigs(reads,ref,stdout=None,stderr=None,threads=1,base="minimap2",workdir="."):
    """Map contig sequences with minimap2.
    Returned paths are absolute"""
    ref = os.path.abspath(ref)
    reads = os.path.abspath(reads)
    with util.chdir(workdir,create=True,to_abs=True) as dh_work:
        _file_name = dh_work.make_file
        bam = _file_name(base, "bam")
        step_dir = "{}.tmp".format(base)
        cmd =  """\
        ln -s {ref} ref.fa
        (minimap2 -x asm10 -a ref.fa {reads} \\
        | samtools view -@ {threads} -bSu - \\
        | samtools sort -@ {threads} -T {reads}.bwa.tmp -o {bam} -) && \\
        samtools index {bam} && \\
        rm -rf {reads}.bwa.tmp
        """
        run_step(cmd,
                 step_dir=step_dir,
                 descr="Map contigs with Minimap2",
                 stdout=stdout, stderr=stderr, args=locals(),
                 shell=True,
                 rm_step_dir=True)
    return dict(bam=bam)


def fasta_index(fasta, stdout=None, stderr=None):
    """Index FASTA file with samtools"""
    cmd = """\
    samtools faidx {}
    """.format(fasta)
    run_step(cmd,
             descr="Index FASTA file with samtools",
             stdout=stdout,
             stderr=stderr)

def picard_dedupe(bam_inp, bam_out, stdout=None, stderr=None, remove_duplicates=False,
                  remove_sequencing_duplicates=True,index_out=True,
                  picard_args="MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=400",
                  metrics_out=None):
    """Mark and optionally remove duplicates with Picard"""
    _picard_bool = lambda x: "true" if  x else "false"
    remove_duplicates = _picard_bool(remove_duplicates)
    remove_sequencing_duplicates = _picard_bool(remove_sequencing_duplicates)
    rm_metrics = False
    if not metrics_out:
        metrics_out = "{}.pic_dd_metrics".format(bam_out)
        rm_metrics= True
    cmd = """\
    picard -XX:ParallelGCThreads=1 MarkDuplicates \
                        {picard_args} \
                        REMOVE_DUPLICATES={remove_duplicates} \
                        REMOVE_SEQUENCING_DUPLICATES={remove_sequencing_duplicates} \
                        INPUT={bam_inp} \
                        OUTPUT={bam_out} \
                        METRICS_FILE={metrics_out}
    """
    run_step(cmd,
             descr="Mark or remove duplicates with Picard MarkDuplicates",
             stdout=stdout, stderr=stderr, args=locals())
    if index_out:
        run_step("samtools index '{}'".format(bam_out),
                 descr="Index BAM post-Picard MarkDuplicates",
                 stdout=stdout, stderr=stderr)
    if rm_metrics:
        os.unlink(metrics_out)
    return dict(bam=bam_out)

def map_and_pilon(inp_reads,contigs,stdout=None,stderr=None,threads=1,gc_threads=None,base="asm",workdir="mpl_asm",
                  bbmap_args="ambig=random minid=0.96 idfilter=0.96 inserttag idtag basecov=basecov.txt",
                  pilon_args="--changes --tracks --iupac --mindepth 5 --fix snps,gaps --iupacminfreq 0.7 --iupacminqualsum 150",
                  do_picard_dedupe=False,
                  picard_dedupe_args="MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=400",
                  remove_duplicates=False,
                  remove_sequencing_duplicates=True,
                  ram=4*1024,
                  rams=2*1024):
    """Map with BBMap and run Pilon.
    Returned paths are absolute"""
    ## Regarding picard deduplication and deduplication in general (e.g. with clumpify):
    ## These are useful posts on dedupping:
    ##Clumpify
    ##https://www.biostars.org/p/225338/
    ##Picard etc
    ##https://wiki.bits.vib.be/index.php/Q%26A_added_during_the_intro_to_NGS_data_analysis#what_are_duplicates_and_what_to_do_about_them_.28an_introduction.21.29
    ##
    ## Clumpify: keep subs at a non-zero value, otherwise will enrich dataset with reads containing errors.
    ##
    ## For a deep viral amplicon sequencing of mixed samples, the links above probably signify that only optical
    ## dedup might be OK, but not PCR dedup, either with clustering (clumpify) or mapping (picard).
    ## Because, if subs > 0, we will be falsely removing as duplicates reads from different molecules
    ## (even both reads 5' coords can be easily identical with a PCR amplicon, and even more when it is processed with Nextera).
    ## If subs=0, we will be enriching for sequencing errors. The assembler will have a harder time
    ## phasing such error-enriched datasets.
    ## Switching do_picard_dedupe=True here did not improve the assembly of 150 PCR amplicon viral samples, and actually
    ## turned one previosly contiguous assembly into a fragmented one.
    ## For Pilon, do not use by default --fix indels - it sometimes breaks assemblies with very large indels

    if not gc_threads:
        gc_threads = threads
    if rams > ram:
        rams = ram
    with util.chdir(workdir,create=True,to_abs=True) as dh_work:
        _read_names = dh_work.make_paired_files_func(base="read",ext="fastq",n=2)
        _file_name = dh_work.make_file
        ## bbmap-generated indexing script adds _sorted suffix:
        bam = _file_name(base + "_sorted", "bam")
        if do_picard_dedupe:
            sam = _file_name(base + "_bbmap", "sam")
            bam_bbmap = _file_name(base + "_bbmap_sorted", "bam")
        else:
            # bbmap file name is the final file name
            sam = _file_name(base, "sam")
            bam_bbmap = bam
        map_dir = "{}_map".format(base)
        ## From: http://bioinformaticsonline.com/blog/view/35621/bbtools-for-bioinformatician
        ## BBMap is always nondeterministic when run in paired-end mode with multiple threads, because
        ## the insert-size average is calculated on a per-thread basis, which affects mapping; and which
        ## reads are assigned to which thread is nondeterministic. The only way to avoid that would be to
        ## restrict it to a single thread (threads=1), or map the reads as single-ended and then fix
        ## pairing afterward (but fix pairing example only pertains to pairing unmapped reads in case of
        ## contaminant filtering, for example). 
        ## Update: now there is a `deterministic` flag (used below).
        ## And we still do not know if ambig=random uses a random or fixed seed.
        bbmap_threads_arg = ""
        if "threads=" not in bbmap_args:
            bbmap_threads_arg = "threads={}".format(threads)
        run_step("bbmap.sh -Xmx{ram}m -Xms{rams}m overwrite=true in='{inp_reads[0]}' in2='{inp_reads[1]}' "
                 "out={sam} ref={contigs} nodisk sampleseed=10 deterministic "
                 "{bbmap_args} {bbmap_threads_arg} bamscript=bbmap_bam.sh && sh bbmap_bam.sh",
                 step_dir=map_dir,
                 descr="Map reads",
                 stdout=stdout, stderr=stderr, args=locals(),
                 shell=True,
                 rm_step_dir=False,
                 env=dict(_JAVA_OPTIONS="-XX:ParallelGCThreads={gc_threads} -XX:+UseParallelOldGC".format(gc_threads=gc_threads)))
        os.unlink(sam)
        if do_picard_dedupe:
            picard_dedupe(bam_bbmap, bam, stdout=stdout, stderr=stderr, remove_duplicates=remove_duplicates,
                          remove_sequencing_duplicates=remove_sequencing_duplicates, index_out=True,
                          picard_args=picard_dedupe_args)
            os.unlink(bam_bbmap)
        pilon_base = "{}.".format(base)
        pilon_dir = "{}_pilon".format(base)
        pilon_out_seq = os.path.join(pilon_dir, pilon_base + ".fasta")
        ## Use 1 thread for determinism - threads is experimental feature in Pilon anyway
        run_step("pilon --threads 1 --genome {contigs} --bam {bam} "
                 "--outdir . --output {pilon_base} {pilon_args}",
                 step_dir=pilon_dir,
                 descr="Use Pilon to polish reference with mapped reads",
                 stdout=stdout, stderr=stderr, args=locals(),
                 ## pilon does not expose ram option, so set it through the _JAVA_OPTIONS
                 env=dict(_JAVA_OPTIONS="-Xmx{ram}m -Xms{rams}m -XX:ParallelGCThreads={gc_threads} -XX:+UseParallelOldGC".format(ram=ram,rams=rams,gc_threads=gc_threads)))
        abspath = os.path.abspath
        return dict(bam=bam,map_dir=abspath(map_dir),pilon_dir=abspath(pilon_dir),
                    pilon_root=abspath(os.path.join(pilon_dir,pilon_base)),
                    pilon_out_seq=abspath(pilon_out_seq))


## import package module and add argh entry points

def _main():
    import argh
    argh.dispatch_commands([
        run_step,
        map_and_pilon,
        bwa_mem_long,
        picard_dedupe,
        sort_reads,
        clean_reads
    ])

if __name__ == "__main__":
    _main()
