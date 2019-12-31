
from __future__ import absolute_import
from __future__ import unicode_literals
from __future__ import print_function
from __future__ import division
from future import standard_library
standard_library.install_aliases()
from builtins import str
from builtins import *
from . import util
from . import workflow_util
from . import yaml_util
from . import jinja_util
from . import num_util
from .pysteps import run_step, map_and_pilon, minimap2_contigs, fasta_index

import os
import shutil
import sys
import tempfile
import csv

def asm_with_ref(
        reads1="read_1.fq",
        reads2="read_2.fq",
        ref="ref.fasta",
        workdir="work",
        threads=1,
        prefix="cluster",
        args=None
    ):
    """Use Pilon to perform "reference-based" assembly.
    This attempts to edit the existing reference from mapping reads.
    NOTE: this is currently incomplete - it will currently leave the
    reference as-is in regions of low coverage. It has to be modified
    by adding a step that maps reads again to the edited sequence (in case
    initial low coverage was due to deletions),
    and then either masks or breaks low coverage regions. There is also some
    option in Pilon to "emit reads for de-novo reassembly" that can be
    explored.
    """

    if(args):
        args = yaml_util.get_arg_as_yaml(args).copy()
    else:
        args = {}

    inp_reads = [ os.path.abspath(_) for _ in [reads1,reads2] ]
    ref = os.path.abspath(ref)

    os.environ["_JAVA_OPTIONS"] = ''

    stdout = None ## goes to default stream

    print("## Starting reference assembler run for prefix {}".format(prefix))

    with util.chdir(workdir,create=True,to_abs=True) as dh_work:
        _read_names = dh_work.make_paired_files_func(base="read",ext="fastq",n=2)
        _file_name = dh_work.make_file
        bbmap_args = "ambig=random inserttag idtag"
        pilon_mindepth = args.get("pilon_mindepth",5)
        pilon_iupacminfreq = args.get("pilon_iupacminfreq", 0.7)
        pilon_iupacminqualsum = args.get("pilon_iupacminqualsum", 150)
        ## do not use --fix indels - it breaks some reasonable assemblies in mixed samples by creating very large deletions (S20 in SA)
        ## do not use --diploid - it supresses ambiguous calls
        pilon_args = ("--changes --tracks --iupac --mindepth {pilon_mindepth} --fix snps,gaps " \
                     "--iupacminfreq {pilon_iupacminfreq} --iupacminqualsum {pilon_iupacminqualsum}").format(**locals())
        base = "{}_asmref".format(prefix)
        mpl_res = map_and_pilon(inp_reads=inp_reads, contigs=ref,
                      stdout=stdout, threads=threads,base=base,workdir="mpl_asmref",
                      bbmap_args=bbmap_args,
                      pilon_args=pilon_args)

    print("## Done reference assembler run for prefix {}".format(prefix))

def post_ariba_asm_with_ref_web(top_dir,
                                target,
                                sample_id,
                                micgent_root="../../../../"):
    env = jinja_util.init_jinja("gene_extractor")
    configFile = "{}_map.json".format(target)
    pageFile = "{}_map.html".format(target)
    #locus = "A01:4000-8000"
    jinja_util.jinja_dump_template(env,"browser_reads.json",
                        out_fn=os.path.join(top_dir,configFile),
                        args=dict(target=target))
    jinja_util.jinja_dump_template(env,"browser_reads.html",
                        out_fn=os.path.join(top_dir,pageFile),
                        args=dict(target=target,
                                  configFile=configFile,
                                  micgent_root=micgent_root,
                                  sample_id=sample_id))

def top_web_pages(top_dir,
                                target="contigs",
                                micgent_root=""):
    env = jinja_util.init_jinja("gene_extractor")
    #locus = "A01:4000-8000"
    jinja_util.jinja_dump_template(env,"browser_contigs.json",
                        out_fn=os.path.join(top_dir,"browser_contigs.json"),
                        args=dict(target=target))
    jinja_util.jinja_dump_template(env,"browser_contigs.html",
                        out_fn=os.path.join(top_dir,"browser_contigs.html"),
                        args=dict(target=target,
                                  configFile="browser_contigs.json",
                                  micgent_root=micgent_root))
    jinja_util.jinja_dump_template(env,"geno_browser_contigs.json",
                        out_fn=os.path.join(top_dir,"geno_browser_contigs.json"))

def inject_iframe_resizer(finp,fout,script_url):
    """Inject IFrameResizer script into the HTML page that will be loaded in an IFrame"""
    from bs4 import BeautifulSoup
    with open(finp, "r",encoding='utf-8') as inp:
        soup = BeautifulSoup(inp,"lxml")
    tag = soup.new_tag("script",type="text/javascript",src="{}".format(script_url))
    soup.body.append(tag)
    with open(fout, "w",encoding='utf-8') as out:
        out.write(soup.prettify(formatter="html"))
    ## this is altenative simplistic code that is probably more efficient for large files,
    ## but it generates an encoding error on multiqc files
    # script_tag = '<script type="text/javascript" src="{}"></script>'.format(script_url)
    # body_close_tag = '</body>'
    # replacement = script_tag+'\n'+body_close_tag
    # if False:
    #     with open(finp,"r") as inp,\
    #         open(fout,"w") as out:
    #         data = inp.read()
    #         out.write(util.rreplace(data,body_close_tag,replacement,1))

def prepare_sor_pack(manifest,seq,
    toc_prov,
    sfx_target,sfx_version,pack_out):
    import pandas as pd
    from Bio import SeqIO
    pack_dir = tempfile.mkdtemp(suffix = 'sor_pack', prefix = 'tmp_', dir = os.getcwd())
    man = workflow_util.load_manifest(manifest)
    assert not man.duplicated("SeqID").any(), "SeqID in the manifest contains duplicated values"
    seq_out = os.path.join(pack_dir,"seq.fasta")
    SeqID_SOR = []
    SeqID_old = []
    with util.open_text_py23(seq,"r") as finp,\
            util.open_text_py23(seq_out,"w") as fout:
        for rec in SeqIO.parse(finp,"fasta"):
            id_old = rec.id
            rec.id = "{}-{}-{}".format(id_old,sfx_target,sfx_version)
            rec.description = ''
            SeqID_SOR.append(rec.id)
            SeqID_old.append(id_old)
            SeqIO.write(rec,fout,"fasta")
    assert (man.shape[0] == len(SeqID_old)) and (man.SeqID == SeqID_old).all(),"Sequence ID mismatch between manifest and FASTA"
    man["SeqID_SOR"] = SeqID_SOR
    manifest_out = os.path.join(pack_dir, "manifest.tsv")
    man.to_csv(manifest_out, index=False, sep="\t")
    ## collect provenance files into a provenance archive inside pack_dir
    prov_dir = os.path.join(pack_dir,"provenance")
    os.makedirs(prov_dir)
    for toc_rec in toc_prov.itertuples():
        shutil.copy(toc_rec.File,os.path.join(prov_dir,toc_rec.File))
    toc_prov.to_csv(os.path.join(prov_dir,"toc.tsv"), index=False, sep="\t")
    util.dir_to_tar(prov_dir, os.path.join(pack_dir,"provenance.tgz"))
    shutil.rmtree(prov_dir, ignore_errors=True)
    util.dir_to_tar(pack_dir, pack_out)
    shutil.rmtree(pack_dir, ignore_errors=True)


def post_extractor(manifest_all="manifest_out_all.tsv",
                   contigs_all="seq_out_all.fasta",
                   manifest="manifest_out.tsv",
                   contigs="seq_out.fasta",
                   manifest_sum="manifest_out_sum.tsv",
                   manifest_dict="manifest_out_dict.tsv",
                   qc_report="multiqc.html",
                   qc_data="multiqc_data.zip",
                   ref_common="ref_common.fasta",
                   micgentjs_tgz="micgentjs.tar",
                   sample_tars="sample_tars.txt",
                   sor_sfx_target="X",
                   sor_sfx_version=1,
                   wf_inputs="wf_inputs.yaml",
                   prepareref_tgz="aribaref.tgz",
                   iframe_resizer_url="https://cdnjs.cloudflare.com/ajax/libs/iframe-resizer/4.2.1/iframeResizer.contentWindow.min.js",
                   out_dir="web",
                   out_tar="web.tar",
                   out_sor="sor_pack.tgz"
                   ):
    import pandas as pd
    pjoin = os.path.join
    ## the command below must create out_dir/micgentjs; such input archive can be created,
    ## for example, if you run this command in a seperate micgentjs Git repo (note the slash in --prefix):
    ## git archive -o micgentjs.tar --prefix=micgentjs/ HEAD && (gzip -c micgentjs.tar > micgentjs.tgz).
    ## In other words, micgentjs directory name should be part of the archive.
    util.tar_to_dir(micgentjs_tgz,out_dir)
    ## therefore, micgent_root is the out_dir itself
    micgent_root = ""

    def _copy_if_diff_name(finp,fout):
        if not (os.path.exists(fout) and os.path.samefile(finp,fout)):
            shutil.copy(finp,fout)

    toc = []
    def _copy_and_describe(finp,fout,ind,descr,out_dir=out_dir,toc=toc,skip_if_missing=True):
        fout_path = pjoin(out_dir, fout)
        if finp and os.path.exists(finp):
            shutil.copy(finp, fout_path)
        if (not skip_if_missing) or os.path.exists(fout_path): 
            toc.append([ind,descr,fout])

    ## collect web report files into the output web dir
    _copy_and_describe(manifest,"manifest_out.tsv","1.1","Final manifest tab-delimited text file")
    _copy_and_describe(contigs,"seq_out.fasta","1.2","Final assemblies FASTA format")
    _copy_and_describe(manifest_all,"manifest_out_all.tsv","2.1","Manifest with all input samples tab-delimited text file")
    _copy_and_describe(contigs_all,"seq_out_all.fasta","2.2","All assemblies FASTA format")
    _copy_and_describe(manifest_sum,"manifest_out_sum.tsv","3.1","Manifest with per-sample summary metrics tab-delimited text file")
    _copy_and_describe(manifest_dict,"manifest_out_dict.tsv","3.2","Manifest data dictionary tab-delimited text file")

    inject_iframe_resizer(qc_report,pjoin(out_dir,"multiqc.html"),script_url=iframe_resizer_url)

    _copy_and_describe(None,"multiqc.html","4.1","MultiQC report HTML format")
    _copy_and_describe(qc_data,"multiqc_data.zip","4.2","MultiQC data ZIP archive")
    _copy_and_describe(ref_common,"ref_common.fasta","5","Common FASTA reference for contig browser")
    _copy_and_describe(wf_inputs,"wf_inputs.yaml","6","Input parameters for the workflow in YAML format")
    _copy_and_describe(prepareref_tgz,"aribaref.tgz","7","Ariba reference pack TAR archive. Contains FASTA reference sequences in the file `inputs.fasta`")
    out_sor_base = "sor_pack.tgz"

    with open(sample_tars, "r") as inp_samp:
        sample_tars =  [ os.path.abspath(fname) for fname in (_.strip() for _ in inp_samp) if fname ]
    print("## Received these sample Web reports: "+str(sample_tars))
    with util.chdir(out_dir, to_abs=True) as dh_work:
        _file_name = dh_work.make_file
        minimap2_contigs("seq_out.fasta", "ref_common.fasta",threads=1, base="seq_out.fasta")
        fasta_index("ref_common.fasta")
        top_web_pages(top_dir=".",micgent_root=micgent_root)
        with util.chdir("samples",create=True):
            for sample_tar in sample_tars:
                sample_id = os.path.basename(sample_tar).rsplit('.asmqc.tar',1)[0]
                util.tar_to_dir(sample_tar,sample_id)
        _copy_and_describe(None, out_sor_base, "8", "TAR archive for loading into System Of Record (SOR)",skip_if_missing=False)
        toc_tab = "toc.tsv"
        with open(toc_tab, "w") as out_toc:
            toc_writer = csv.writer(out_toc, dialect='excel-tab', lineterminator='\n')
            toc_writer.writerow(["Index","Description","File"])
            for rec in toc:
                toc_writer.writerow(rec)

        toc = pd.read_csv(toc_tab,dialect="excel-tab")
        with open("toc.html","w",encoding="utf-8") as out:
            out.write(toc.to_html(escape=False,index=False,formatters=dict(File=lambda x: '<a href="{}">{}</a>'.format(x,x))))
        toc_prov = toc[~toc.Description.isin(("TAR archive for loading into System Of Record (SOR)","Final manifest","Final assemblies"))]
        prepare_sor_pack(manifest=manifest,seq=contigs,
                           toc_prov=toc_prov,
                           sfx_target=sor_sfx_target,sfx_version=sor_sfx_version,
                           pack_out=out_sor_base)

    _copy_if_diff_name(os.path.join(out_dir,out_sor_base),out_sor)
    util.dir_to_tar(out_dir,out_tar)


def filter_assemblies(manifest="manifest.tsv",
                   contigs="seq.fasta",
                   manifest_out="manifest_out.tsv",
                   contigs_out="seq_out.fasta",
                   manifest_out_all="manifest_out_all.tsv",
                   manifest_out_sum="manifest_out_sum.tsv",
                   manifest_out_dict="manifest_out_dict.tsv",                   
                   args=None
                   ):
    """Filter the combined set of assemblies after the initial extraction workflow
    """
    import pandas as pd
    import numpy as np
    from Bio import SeqIO

    if (args):
        args = yaml_util.get_arg_as_yaml(args).copy()
    else:
        args = {}
    man = workflow_util.load_manifest(manifest)
    policy = args.get("policy",None)
    man_flt = man[man.SeqID.notna()].copy()
    man_flt.set_index("SeqID", drop=False, append=False, inplace=True, verify_integrity=True)
    ## add various columns that may be computed by any filtering policies, so that our output manifest header is stable
    flt_metrics_cols = dict(ctg_cov_ratio=np.nan,ctg_ref_bases_ratio=np.nan,SeqCountPostFilter=0,DroppedMultiseq=False)
    man_flt = num_util.create_columns_in_df(man_flt,flt_metrics_cols)
    man = num_util.create_columns_in_df(man,flt_metrics_cols)

    ## this will accumulate various computed filtering metrics, saved here as {metric name => DataFrame([SeqID,metric name])}
    ## before the corresponding filter is applied (since it can drop records, but we want to keep the metrics)
    flt_metrics = {}
    if man_flt.shape[0] > 0:
        if policy:
            if policy == "default":
                ctg_ref_bases_min = args.get("ctg_ref_bases_min",0)
                ctg_ref_bases_ratio_min = args.get("ctg_ref_bases_ratio_min",0)
                ctg_len_min = args.get("ctg_len_min",0)
                man_flt["ctg_ref_bases_ratio"] = man_flt.Asm_ref_base_assembled/man_flt.Asm_ref_len
                flt_metrics["ctg_ref_bases_ratio"] = man_flt[["SeqID","ctg_ref_bases_ratio"]].set_index("SeqID",verify_integrity =True)
                man_flt = man_flt[(man_flt.Asm_seq_len >= ctg_len_min) & \
                    (man_flt.Asm_ref_base_assembled >= ctg_ref_bases_min) & \
                    (man_flt.ctg_ref_bases_ratio >= ctg_ref_bases_ratio_min)].copy()
                if man_flt.shape[0] > 0:
                    ## with no rows, the next line raises
                    ctg_cov_ratio = man_flt.groupby("SampleID").Asm_ctg_cov.transform(lambda x: x / x.sum())
                    assert (man_flt.index == ctg_cov_ratio.index).all(),"Detected apparent bug in Pandas - order lost after transform()"
                    man_flt["ctg_cov_ratio"] = ctg_cov_ratio
                    flt_metrics["ctg_cov_ratio"] = man_flt[["SeqID","ctg_cov_ratio"]].set_index("SeqID",verify_integrity =True)
                    ctg_cov_ratio_min = args.get("ctg_cov_ratio_min", 0)
                    ctg_cov_min = args.get("ctg_cov_min", 0)
                    man_flt = man_flt[(man_flt.ctg_cov_ratio >= ctg_cov_ratio_min) & \
                        (man_flt.Asm_ctg_cov >= ctg_cov_min)]
                if man_flt.shape[0] > 0:
                    ## the two lines below seems like the least verbose way of assigning group count to each original ungrouped row
                    man_flt["SeqCountPostFilter"] = 1
                    seq_count = man_flt.groupby("SampleID").SeqCountPostFilter.transform(lambda x: x.sum())
                    assert (man_flt.index == seq_count.index).all(),"Detected apparent bug in Pandas - order lost after transform()"
                    ## count after all filters but before rejecting mixed samples
                    man_flt["SeqCountPostFilter"] = seq_count
                    flt_metrics["SeqCountPostFilter"] = man_flt[["SampleID","SeqCountPostFilter"]].drop_duplicates().set_index("SampleID",verify_integrity =True)
                    multiseq_policy = args.get("multiseq_policy","accept")
                    if multiseq_policy == "accept":
                        pass
                    elif multiseq_policy == "reject":
                        man_flt["DroppedMultiseq"] = man_flt.SeqCountPostFilter>1
                        flt_metrics["DroppedMultiseq"] = man_flt[["SeqID","DroppedMultiseq"]].set_index("SeqID",verify_integrity =True)
                        man_flt = man_flt[~man_flt.DroppedMultiseq]
                    else:
                        raise ValueError("Unknown value of multiseq_policy: {}".format(multiseq_policy))
            else:
                raise ValueError("Unknown filtering policy name: '{}'".format(policy))
    seq_id_keep = set(man_flt.index)
    seq_id_out = []
    with util.open_text_py23(contigs,"r") as finp,\
            util.open_text_py23(contigs_out,"w") as fout:
        for rec in SeqIO.parse(finp,"fasta"):
            if rec.id in seq_id_keep:
                SeqIO.write(rec,fout,"fasta")
                seq_id_out.append(rec.id)
    assert (man_flt.shape[0] == len(seq_id_out)) and (man_flt.SeqID == seq_id_out).all(),"Sequence ID mismatch between output manifest and FASTA"
    man_flt.drop(columns=["DroppedMultiseq"],inplace=True)
    man_flt = num_util.add_column_name_prefix(man_flt,"Asm_",flt_metrics_cols.keys())
    man_flt.to_csv(manifest_out, index=False, sep="\t")
    man_bak = man.copy()
    ## make sure that the manifest has all metric fields - update with the accumulated metric df if exists,
    ## after assigning the default value
    for flt_c in flt_metrics_cols.keys():
        if flt_c in flt_metrics:
            man.set_index(flt_metrics[flt_c].index.name,drop=False,inplace=True)
            man.update(flt_metrics[flt_c])
    
    assert num_util.keys_equal(man,man_bak,["SeqID"]), "Missing or extra rows after join"

    man = man.reset_index(drop=True)
    man = num_util.add_column_name_prefix(man,"Asm_",flt_metrics_cols)

    man["Asm_Failed"] = ~man.SampleID.isin(man_flt.SampleID)
    man["Asm_Dropped"] = (~man.SeqID.isin(man_flt.SeqID)) & man.SeqID.notnull()
    man.to_csv(manifest_out_all, index=False, sep="\t")
    ## could not shake the float type that appeared for unknown reason, using
    ## explicit conversion
    man["Asm_Dropped"] = man.Asm_Dropped.astype(int)

    ## A section that computes a sample-level summary status table
    man_gb = man.groupby("SampleID")
    ## count() only counts non-null values
    man_sm = man_gb.SeqID.count().to_frame(name='Asm_SeqCountPreFilter') #converts series to one-column df with index SampleID

    ## join with each aggregated metric separately to avoid NA in some metric dropping whole rows
    ## and affecting other metrics
    man_sm = (man_sm
        .join(man_gb.agg({'Asm_Failed':np.all}))
        .join(man_gb.agg({'Asm_Dropped':"sum"})).rename(columns={'Asm_Dropped':'Asm_DroppedCount'})
        .join(man_gb.agg({'Asm_DroppedMultiseq':np.any}))
        .join(man_gb.agg({'Asm_SeqCountPostFilter': "max"})) # max because some values can be 0 if SeqID was filtered earlier
    ).reset_index() #reset_index converts SampleID index to column
    man_sm["Asm_SeqCountPostFilter"] = man_sm.Asm_SeqCountPostFilter.astype(int)
    man_sm["Asm_SeqCountFinal"] = man_sm.Asm_SeqCountPreFilter - man_sm.Asm_DroppedCount
    #import pdb; pdb.set_trace()
    assert num_util.keys_equal(man_sm[man_sm.Asm_SeqCountFinal>0],
        man_flt.groupby("SampleID").SeqID.count().to_frame("Asm_SeqCountFinal").reset_index(),
        keys=["SampleID","Asm_SeqCountFinal"]),"Sequence counts computed from the sample summary table and the output manifest do not match"
    man_sm.to_csv(manifest_out_sum, index=False, sep="\t")
    man_dict = [
    ("SampleID","ID of the sequencing sample (unique record ID in the input manifest and in the output summary 'manifest_out_sum.tsv')"),
    ("SeqID","ID of the output sequence (unique record ID in the output 'manifest_out.tsv' and 'manifest_out_all.tsv')"),
    ("SeqID_SOR","Sequence ID that exactly matches FASTA ID in the System Of Record (SOR) pack (with SOR FASTA suffix)"),    
    ("file1","Absolute file name of the first read in a pair"),
    ("file2","Absolute file name of the second read in a pair"),
    ("Asm_cluster","ID of the read cluster forming the assembly unit"),    
    ("Asm_ref_name","Sequence ID of the best matching recruiting reference for a given cluster"),
    ("Asm_reads","Number of reads in the assembly cluster"),
    ("Asm_ref_len","Length of the reference"),
    ("Asm_ref_base_assembled","Number of bases in the reference that are present in the assembly"),
    ("Asm_pc_ident","Percentage sequence identity between assembly and reference"),
    ("Asm_ctg_len","Length of the assembled contig"),
    ("Asm_ctg_cov","Median read coverage depth of the assembled contig"),
    ("Asm_seq_len","Length of the final sequence extracted from the assembled contig"),
    ("Asm_SeqStatus","Type of output sequence extracted: [gene - ORF matching reference gene; match - output trimmed to matching reference, but is not a gene; asm - whole contig"),
    ("Asm_msg","'Finished' if no warnings were returned by Ariba assembly, otherwise concatenated error or warning messages. An error for one cluster can be a normal status if other cluster has assembled OK."),
    ("Asm_sig_inp1","Blake2b signature of the first read file in a pair, computed with a key 'sig_inp_key' in 'wf_inputs.yaml' file."),
    ("Asm_sig_inp2","Blake2b signature of the second read file in a pair, computed with a key 'sig_inp_key' in 'wf_inputs.yaml' file."),    
    ("Asm_ctg_cov_ratio","Ratio of median read coverage depth of the assembled contig to the sum of such depths of all contigs"),
    ("Asm_ctg_ref_bases_ratio","Ratio of reference bases covered by the best matching contig in assembly"),
    ("Asm_SeqCountFinal","Number of sequences in the final manifest, per sample ID"),
    ("Asm_SeqCountPreFilter","Number of sequences per sample that existed before the initial filter"),    
    ("Asm_SeqCountPostFilter","Number of sequences per sample that passed the initial filter"),
    ("Asm_DroppedMultiseq","Flag any sequences that passed the initial filter but were dropped due to policy of handling samples with multiple outputs"),
    ("Asm_Failed","Flag samples that did not result in any final sequences, either due to no sequences assembled or all sequences filtered out"),
    ("Asm_Dropped","Flag sequences that were present in the assemblies but did not pass to the final manifest due to filtering"),
    ("Asm_DroppedCount","Number of sequences per sample with a flag Asm_Dropped")
    ]
    man_dict = pd.DataFrame.from_records(man_dict,columns=["Variable","Description"])
    man_dict.to_csv(manifest_out_dict, index=False, sep="\t")


def post_ariba_asm_with_ref(
        ariba_dir="ariba",
        out_dir="post_asm",
        out_tar="post_asm.tar",
        threads=1,
        sample_id="UnspecifiedSample",
        args=None,
        leave_out_dir=False,
        basecov_asm=None,
        basecov_ref=None
    ):
    """Use Pilon to perform "reference-based" assembly from Ariba run directory.
    Edit de-novo assembled contigs and generate diagnostic files.
    """
    import pandas as pd

    if(args):
        args = yaml_util.get_arg_as_yaml(args).copy()
    else:
        args = {}

    if not basecov_asm:
        basecov_asm = os.path.abspath("{}.basecov_asm.txt".format(sample_id))
    if not basecov_ref:
        basecov_ref = os.path.abspath("{}.basecov_ref.txt".format(sample_id))

    os.environ["_JAVA_OPTIONS"] = '-XX:ParallelGCThreads=1'

    stdout = sys.stdout
    stderr = sys.stderr
    ariba_dir = os.path.abspath(ariba_dir)
    pjoin = os.path.join
    report_file = pjoin(ariba_dir,"report.tsv")
    rep = pd.read_table(report_file)

    #rep_aggr = rep.groupby(["ref_name", "ctg"], as_index=False).head(1)
    clusters = rep.cluster.unique()
    out_dir_clusters = pjoin(out_dir,"clusters")
    os.makedirs(out_dir_clusters)
    basecov_asm_df = []
    basecov_ref_df = []
    for cluster in clusters:
        top_workdir = pjoin(out_dir_clusters,cluster)
        cluster_dir = pjoin(ariba_dir,"clusters",cluster)
        inp_reads = [pjoin(cluster_dir,"reads_{}.fq".format(i_read)) for i_read in (1,2)]
        contigs_inp = pjoin(cluster_dir,"assembly.fa")
        ref_inp = pjoin(cluster_dir,"reference.fa")
        with util.chdir(top_workdir, create=True, to_abs=True) as dh_work:
            _file_name = dh_work.make_file
            contigs = _file_name("assembly","fa")
            ref = _file_name("reference", "fa")
            shutil.copy(contigs_inp,contigs)
            shutil.copy(ref_inp, ref)

            run_step("samtools faidx {}".format(ref)) #for genome browser
            run_step("samtools faidx {}".format(contigs)) #for genome browser
            minimap2_contigs(reads=ref, ref=contigs, stdout=stdout, threads=threads, base="ref_to_asm", workdir=".")
            minimap2_contigs(reads=contigs, ref=ref, stdout=stdout, threads=threads, base="asm_to_ref", workdir=".")

            workdir = "mpl_asm"
            print("## Mapping and running Pilon on the assembly for cluster {}".format(cluster))
            ## bbmap often ignores threads=8 argument, uses 46 cores and takes extremely long to run;
            ## trying to restrict to a single thread here since on this short reference
            ## it should be fast enough even with deep coverage
            bbmap_args = "ambig=random minid=0.96 idfilter=0.96 inserttag idtag basecov=basecov.txt"
            #do not use --fix indels - it sometimes breaks assemblies with very large indels
            pilon_args = "--changes --tracks --iupac --mindepth 5 --fix snps,gaps --iupacminfreq 0.7 --iupacminqualsum 150"
            base = "asm"
            mpl_res = map_and_pilon(inp_reads=inp_reads, contigs=contigs,
                          stdout=stdout,stderr=stderr,threads=threads,base=base,workdir=workdir,
                          bbmap_args=bbmap_args,
                          pilon_args=pilon_args)
            basecov = os.path.join(mpl_res["map_dir"],"basecov.txt")
            if os.path.exists(basecov):
                basecov_asm_df.append(pd.read_table(basecov))

            workdir = "mpl_ref"
            print("## Mapping and running Pilon on the reference for cluster {}".format(cluster))
            bbmap_args = "ambig=random inserttag idtag local basecov=basecov.txt" ## `local` needed when the reference is longer than the amplicon
            ## we do not allow Pilon to break the reference (--fix break) because we observe it breaking it at coverage changes in
            ## amplicon datasets; same happens with --fix indels (deletes large chunks of sequence)
            pilon_args = "--changes --tracks --iupac --mindepth 5 --fix snps,gaps"
            base = "asm"
            mpl_res = map_and_pilon(inp_reads=inp_reads, contigs=ref,
                          stdout=stdout,stderr=stderr,threads=threads,base=base,workdir=workdir,
                          bbmap_args=bbmap_args,
                          pilon_args=pilon_args)
            basecov = os.path.join(mpl_res["map_dir"],"basecov.txt")
            if os.path.exists(basecov):
                basecov_ref_df.append(pd.read_table(basecov))

            minimap2_contigs(reads=mpl_res["pilon_out_seq"], ref=ref, stdout=stdout, threads=threads, base="asm_ref_to_ref", workdir=".")

            print("## Done mapping to contigs and to the reference assembler for cluster {}".format(cluster))
            for target in ("ref","asm"):
                post_ariba_asm_with_ref_web(top_dir="",target=target,sample_id=sample_id)
            print("## Generated Web pages for the reference assembler for cluster {}".format(cluster))

    if len(basecov_asm_df) > 0:
        basecov_asm_df = pd.concat(basecov_asm_df,axis=0)
    else:
        basecov_asm_df = pd.DataFrame(columns=["#RefName","Pos","Coverage"])
    basecov_asm_df.to_csv(basecov_asm,index=False, sep="\t")
    if len(basecov_ref_df) > 0:
        basecov_ref_df = pd.concat(basecov_ref_df,axis=0)
    else:
        basecov_ref_df = pd.DataFrame(columns=["#RefName","Pos","Coverage"])
    basecov_ref_df.to_csv(basecov_ref,index=False, sep="\t")

    util.dir_to_tar(out_dir,out_tar)
    print("## Created a tar archive with Web pages: {}".format(out_tar))
    if not leave_out_dir:
        util.rmtree(out_dir,ignore_errors=True)



## import package module and add argh entry points

def _main():
    import argh
    argh.dispatch_commands([
        post_ariba_asm_with_ref,
        post_extractor,
        filter_assemblies
    ])

if __name__ == "__main__":
    _main()
