from __future__ import absolute_import
from __future__ import unicode_literals
from __future__ import print_function
from __future__ import division
from future import standard_library
standard_library.install_aliases()
from builtins import str
from builtins import *
from . import util
from . import yaml_util
from .pysteps import run_step, replace_fasta_ids, map_and_pilon, sort_reads

from subprocess import check_call
import os
import shlex
import shutil
import sys
import itertools
import traceback
import subprocess

def asm_for_skewed_coverage(
        reads1="read_1.fq",
        reads2="read_2.fq",
        contigs="contigs.fasta",
        workdir="work",
        threads=1,
        log="asm.log",
        prefix="cluster",
        extra_args=None
    ):
    """De-novo script that conforms to the `plugin` assembler interface for Ariba.
    This particular script targets datasets with highly skewed uneven coverage.
    It performs trimming to a set of lengths, followed digital depth normalization, 
    followed by Spades assembly, followed by Pilon polishing. 
    This is done for a grid of length and depth values.

    :param extra_args: YAML string that encode a map with this optional keys:
        - deterministic: If true, make extra precautions to ensure deterministic
          behaviour
        - n_iter_retry: Try each assembly iteration up to this number of times in order
          to avoid non-determinism due to random external sources of those 
          failures that we supress in the exception handler
          Note: this will apply even if `deterministic` is False. If `deterministic` is
          True, this will be automatically set to be no less than 3

    """
    ## Note - it is not possible to use multiprocessing.Pool here because this is
    ## already in a worker process of a Pool, but it is possible to use
    ## ThreadPool, which might be even more efficient if the threads mostly call
    ## executables in the subprocesses. See discussion here:
    ## https://stackoverflow.com/questions/6974695/python-process-pool-non-daemonic/34069003

    if(extra_args):
        args = yaml_util.get_arg_as_yaml(extra_args).copy()
    else:
        args = {}

    deterministic = args.get("deterministic",False)

    n_iter_retry = args.get("n_iter_retry",1)
    if deterministic:
        n_iter_retry = max(n_iter_retry,3)

    inp_reads = [ os.path.abspath(_) for _ in [reads1,reads2] ]

    os.environ["_JAVA_OPTIONS"] = '-XX:ParallelGCThreads=1 -XX:+UseParallelOldGC' #'-XX:-UseCompressedOops'

    workdir = os.path.abspath(workdir)
    os.makedirs(workdir,exist_ok=True)

    with open(log, "w") as stdout:
        stderr = stdout
        with util.chdir(workdir,create=True,to_abs=True) as dh_work:
            _read_names = dh_work.make_paired_files_func(base="read",ext="fastq",n=2)
            _file_name = dh_work.make_file

            print("## Running on the node={} with environment={}".format(os.uname(),sorted(os.environ.items())),file=stdout)

            out_seqs = []
            norm_target_list = args.get("norm_target_list",(50,100,200))
            minlen_list = args.get("minlen_list",(50,128,200))
            #norm_target_list = (50,100)
            #minlen_list = (128,)

            for norm_target, minlen in itertools.product(norm_target_list,minlen_list):
                workdir_iter = os.path.join(workdir,"iter_{}_{}".format(norm_target,minlen))
                for i_try in range(n_iter_retry):
                    try:
                        out_seqs.append(asm_iter(inp_reads=inp_reads, workdir=workdir_iter,
                                                threads=threads, stdout=stdout, stderr=stderr, args=args,
                                                norm_target=norm_target,minlen=minlen,
                                                prefix=prefix,
                                                deterministic=deterministic))
                    except:
                        if i_try < (n_iter_retry-1):
                            msg = " Retrying..."
                        else:
                            msg = " Moving on..."

                        print("## Ignoring exception raised during assembly iteration for norm_target={} and minlen={}.{}".\
                              format(norm_target,minlen,msg),file=stdout)
                        traceback.print_exc(file=stdout)
                        try:
                            print(subprocess.run("vmstat",stdout=subprocess.PIPE,encoding="utf8").stdout,file=stdout)
                        except:
                            pass
                    else:
                        break

            dedupe_seq = _file_name("dedupe","fasta")
            ## Ariba will pick the best assembly based on coverage of reference,
            ## dedupe can knock off shorter sequences from the overall "best"
            ## assembly, so skip dedupe for now and simply concatenate the assemblies
            ## into one file
            if False:
                out_seqs_comma = ','.join(out_seqs)
                run_step("dedupe.sh -Xmx8g -Xms2g overwrite=true in={out_seqs_comma} "
                         "out={dedupe_seq} threads={threads} "
                         "minidentity=99",
                         step_dir="asm_dedupe",
                         descr="Dedupe assemblies from different iterations",
                         stdout=stdout, stderr=stderr, args=locals())
            else:
                if out_seqs:
                    check_call(["sh", "-c", "cat {} > {}".format(' '.join(out_seqs), dedupe_seq)])
                else:
                    with open(dedupe_seq,"w"):
                        pass

        shutil.copyfile(dedupe_seq, contigs)

        print("## Done plugin assembler run",file=stdout)


def asm_iter(inp_reads,workdir,threads,stdout,stderr,args,norm_target,minlen,prefix,deterministic):

    workdir = os.path.abspath(workdir)

    print("## Starting assembly iteration for norm_target={} and minlen={}".format(norm_target,minlen),file=stdout)

    with util.chdir(workdir,create=True,to_abs=True) as dh_work:
        _read_names = dh_work.make_paired_files_func(base="read",ext="fastq",n=2)
        _file_name = dh_work.make_file

        def sort_reads_cond(inp_reads, out_sfx, base, 
            threads=threads,stdout=stdout,stderr=stderr, deterministic=deterministic):
            if deterministic:
                out_reads = _read_names(out_sfx)
                sort_reads(inp_reads,out_reads,threads=threads,stdout=stdout,stderr=stderr,base=base)
            else:
                out_reads = inp_reads
            return out_reads

        trim_uns_reads = _read_names("tr_uns")
        run_step("bbduk.sh -eoom -Xmx4g  -Xms2g overwrite=true in='{inp_reads[0]}' in2='{inp_reads[1]}' "
                 "out={trim_uns_reads[0]} out2={trim_uns_reads[1]} threads={threads} "
                 "minlen={minlen} ordered",
                 step_dir="trim",
                 descr="Quality-trim and length-filter reads",
                 stdout=stdout, stderr=stderr, args=locals())

        trim_reads = sort_reads_cond(trim_uns_reads,"tr","trim_srt")

        read_correction = args.get("read_correction","ext")
        if not read_correction in ("ext","spades","none"):
            raise ValueError("Unknown value of read_correction: {}".format(read_correction))

        if read_correction == "ext":
            ecorr_uns_reads = _read_names("ec_uns")
            run_step("tadpole.sh -eoom -Xmx4g -Xms2g overwrite=true in={trim_reads[0]} in2={trim_reads[1]} "
                     "out={ecorr_uns_reads[0]} out2={ecorr_uns_reads[1]} threads={threads} "
                     "mode=correct k=62 prefilter=2 prepasses=auto",
                     step_dir="ecorr",stderr=stderr, 
                     descr="Error-correct reads",
                     stdout=stdout, args=locals())
            ecorr_reads = sort_reads_cond(ecorr_uns_reads,"ec","ec_srt")

        else:
            ecorr_reads = trim_reads
        
        if norm_target > 0:

            norm_uns_reads = _read_names("nr_uns")
            ## using prefilter=t is crtitical for memory and speed:
            ## despite what the manual says, v.37.99 gives out-of-heap exception
            ## if given less than 24G heap on typical RSV amplicon data, and runs
            ## for a long time even with 24G. With prefilter, it fits into any ram
            ## (as promised) and is much faster
            run_step("bbnorm.sh -eoom -Xmx4g -Xms2g overwrite=true in='{ecorr_reads[0]}' in2='{ecorr_reads[1]}' "
                     "out={norm_uns_reads[0]} out2={norm_uns_reads[1]} threads={threads} "
                     "target={norm_target} prefilter mindepth=2 passes=2 fixspikes deterministic",
                     step_dir="norm",
                     descr="Digitally normalize depth of coverage",
                     stdout=stdout, args=locals())
            norm_reads = sort_reads_cond(norm_uns_reads,"nr","norm_srt")
        else:
            norm_reads = ecorr_reads

        unmerged_uns_reads = _read_names("um_uns")
        merged_uns_reads = _read_names("m_uns", n=1)
        run_step("bbmerge.sh -eoom -Xmx4g -Xms2g overwrite=true in='{norm_reads[0]}' in2='{norm_reads[1]}' "
                 "outu1={unmerged_uns_reads[0]} outu2={unmerged_uns_reads[1]} threads={threads} "
                 "out={merged_uns_reads} vstrict ordered",
                 step_dir="merge",
                 descr="Merge overlapping paired end reads",
                 stdout=stdout, stderr=stderr, args=locals())

        if deterministic:
            unmerged_reads = _read_names("um")
            merged_reads = _read_names("m", n=1)
            sort_reads([merged_uns_reads], [merged_reads], threads=threads, stdout=stdout, stderr=stderr, base="merge_srt")
            sort_reads(unmerged_uns_reads, unmerged_reads, threads=threads, stdout=stdout, stderr=stderr, base="merge_un_srt")
        else:
            unmerged_reads = unmerged_uns_reads
            merged_reads = merged_uns_reads
        
        spades_mode = args.get("spades_mode", "sc")
        spades_out_seq_base = "contigs.fasta"
        spades_opt = ""
        if spades_mode == "rna":
            spades_opt = "--rna -k 33,55,77,99,127"
            spades_out_seq_base = "transcripts.fasta"
        elif spades_mode in ("sc", "wgs"):
            ##still need to use Hammer read correction in Spades (not --only-assembler), because prior to that we do it
            ##with Tadpole on reads selected by length but before diginorm, in order to speed up bbnorm
            ## Somehow turning it off here results in a couple of fragmented assemblies.
            ## Keeping Spades internal polishing for now (--careful).
            ##It would be nice disabling RR in view of criticism in Unicycler paper, and hoping that Pilon will do some RR downstream,
            ##but this leads to a few fragmented assemblies. Probably, RR is essential component of Spades pipeline in
            ##the SC mode (see their paper on the RR).
            spades_opt = "-k 33,55,77,99,127 --careful"
            #spades_opt = "-k 33,55,77,99,127 --only-assembler --careful"
            #spades_opt = "-k 127 --only-assembler --disable-rr"
      
        if spades_mode == "sc":
            spades_opt = "--sc " + spades_opt
        spades_opt += " --cov-cutoff off"
        spades_dir = "spades"
        spades_out_dir_base = "out"
        spades_out_dir = os.path.join(spades_dir, spades_out_dir_base)
        spades_out_seq = os.path.abspath(os.path.join(spades_out_dir, spades_out_seq_base))
        spades_threads = threads
        if deterministic:
          spades_threads = 1
        ## Spades keeps spending a long time with low cpu utilization when ran with mutiple threads.
        ## Try single-threaded always if that persists?
        ## Spades appears to be deterministic with a fixed number of threads, based on testing and
        ## as confirmed by the developers https://github.com/ablab/spades/issues/111.
        ## Update: I actually caught on non-determinism on very low coverage (~1x)
        ## single-contig assemblies. Every ~1/100 replicates, it reverse-complements
        ## the contig, and also reports slightly different coverage in its name.
        ##
        ## The --pe1 --s2 combination below is designed to provide unmerged and
        ## merged reads as separate libraries as per:
        ## https://twitter.com/spadesassembler/status/907714056387252225
        ## See also https://github.com/tseemann/shovill/blob/master/bin/shovill
        ## TODO: port k-mer range estimation from shovill
        run_step("spades.py -t {spades_threads}  {spades_opt} --pe1-1 {unmerged_reads[0]} --pe1-2 {unmerged_reads[1]} "
                 "--s2 {merged_reads} -o {spades_out_dir_base}",
                 step_dir=spades_dir,
                 descr="Assemble with Spades",
                 stdout=stdout, stderr=stderr, args=locals())
        ## Note: deterministic is always added by pystep
        bbmap_args = "ambig=random minid=0.96 idfilter=0.96 inserttag idtag"

        pilon_mindepth = args.get("pilon_mindepth",5)
        pilon_iupacminfreq = args.get("pilon_iupacminfreq", 0.7)
        pilon_iupacminqualsum = args.get("pilon_iupacminqualsum", 150)
        ## do not use --fix indels - it breaks some reasonable assemblies in mixed samples by creating very large deletions (S20 in SA)
        ## do not use --diploid - it supresses ambiguous calls
        pilon_args = ("--changes --tracks --iupac --mindepth {pilon_mindepth} --fix snps,gaps " \
                     "--iupacminfreq {pilon_iupacminfreq} --iupacminqualsum {pilon_iupacminqualsum}").format(**locals())

        mpl_res = map_and_pilon(inp_reads=inp_reads, contigs=spades_out_seq,
                      stdout=stdout, stderr=stderr, threads=threads, base="asm", workdir="mpl_asm",
                      bbmap_args=bbmap_args,
                      pilon_args=pilon_args)

        fasta_out = os.path.abspath("contigs.fasta")
        replace_fasta_ids(mpl_res["pilon_out_seq"], fasta_out, prefix=prefix, key1=norm_target, key2=minlen)

    print("## Done assembly iteration for norm_target={} and minlen={}".format(norm_target,minlen),file=stdout)
    return fasta_out


## import package module and add argh entry points

def _main():
    import argh
    argh.dispatch_commands([
        asm_for_skewed_coverage
    ])

if __name__ == "__main__":
    _main()
