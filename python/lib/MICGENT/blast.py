"""NCBI BLAST runners and parsers"""
from __future__ import absolute_import
from __future__ import division
from __future__ import unicode_literals
from __future__ import print_function

from future import standard_library
standard_library.install_aliases()
from builtins import str
from builtins import zip
from builtins import *
from past.utils import old_div
from . import util

from MGT.Logging import *
from subprocess import check_call
import csv
import glob, os


def makeblastdb(inp=None,inp_list=None,out=None,dbtype="nucl",**kw):
    """Make BLAST DB.
    @param inp can be a glob or list of globs
    """
    args = []
    for arg_pair in [ ["-"+item[0],item[1]] for item in list(kw.items()) ]:
        args += arg_pair
    inp = list(util.glob_files(files_globs=inp,files_globs_list=inp_list))
    assert len(inp), "No input files provided to makeblastdb"
    inp = " ".join(inp)
    cmd = ["makeblastdb", "-dbtype", dbtype, "-in", inp] + \
          (["-out", out] if out else []) + \
                                                     args
    log.debug(str(cmd))
    check_call(cmd)

def get_blast_dbtype_from_app_name(name):
    ret = "prot"
    if name.endswith("n"):
        ret = "nucl"
    return ret

def get_blast_app(name="blastn"):
    import Bio.Blast.Applications
    return (getattr(Bio.Blast.Applications,"Ncbi{}Commandline".format(name)))


def blast(seq_query=None,
          seq_query_list=None,
          seq_db=None,
          seq_db_list=None,
          app="blastn",
          db_format="fasta",
          blast_out=None,
          dbtype=None,
          makedb_par={},
          blast_par={"outfmt":5},
          short_reference=False,
          return_records=False,
          cleanup=True):
    """Run BLAST, optinally creating a DB.
    Both seq_query and seq_db can be lists of shell globs
    (unless seq_db is already a BLAST DB, in which case it should
    be a single DB alias)
    """
    is_tmp_db = False
    is_tmp_blast_out = False

    if db_format != "blastdb":
        blastdb = util.make_work_file("blastdb")
        if dbtype is None:
            dbtype = get_blast_dbtype_from_app_name(app)
        makeblastdb(inp=seq_db,inp_list=seq_db_list,
                    out=blastdb,dbtype=dbtype,**makedb_par)
        is_tmp_db = True
    else:
        blastdb = seq_db

    if not blast_out:
        blast_out = util.make_work_file("blast_out")
        is_tmp_blast_out = True

    seq_query = list(util.glob_files(files_globs=seq_query,
                                     files_globs_list=seq_query_list))
    if len(seq_query) > 1:
        seq_query_inp = util.make_work_file("blastquery")
        ## blast executables can read query as input stream,
        ## so in principle we can stream from shell `cat`, but
        ## Biopython wrapper needs a real file.
        util.cat_files(seq_query,seq_query_inp)
    else:
        seq_query_inp = seq_query[0]
    blastapp = get_blast_app(app)
    blast_par = blast_par.copy()
    xml_fmt = 5
    blast_par.setdefault("outfmt",xml_fmt)
    if app == "blastn":
        blast_par.setdefault("task", "blastn") #megablast used to be a deafult task in blast+
        if short_reference:
            ##probably, task="blastn-short" would be enough, but for good measure we are adding
            ##specific settings from here: https://www.biostars.org/p/47203/#47207
            blast_par.update(dict(task="blastn-short",word_size=7,dust="no",evalue=10000))
    blast_cln = blastapp(query=seq_query_inp, db=blastdb, out=blast_out, **blast_par)
    log.info(str(blast_cln))
    stdout, stderr = blast_cln()
    if stderr:
        log.warn(stderr)
    record_stream = None
    records = None
    if return_records:
        if blast_par["outfmt"] != xml_fmt:
            raise ValueError("I can only return records for xml BLAST output")
        from Bio.Blast import NCBIXML
        record_stream = open(blast_out,"r")
        records = NCBIXML.parse(record_stream)

    if cleanup:
        if is_tmp_db:
            os.remove(blastdb) #always empty file
            for fname in glob.glob(blastdb+".*"):
                os.remove(fname)
        if is_tmp_blast_out and return_records:
            ## cannot unlink open file on Windows
            try:
                os.remove(blast_out)
                blast_out = None
            except:
                pass

    return dict(stream=record_stream,
                records=records,
                filename=blast_out)

def hsp_filter_default(rec=None,alignment=None,hsp=None):
    """Default filter and converter of BLAST HSP records into tabular (iterator of tuples) structure."""
    if rec is None:
        return ("sbjct_id", "query_id", "pct_ident", "query_cov",
                "query_length", "sbjct_length", "identities", "align_length",
                'query_start','query_end','sbjct_start','sbjct_end',
                'expect','query_strand','sbjct_strand')
    nongap_align_length = hsp.align_length - hsp.gaps
    query_cov = (old_div(float(nongap_align_length), rec.query_length)) if rec.query_length > 0 else 0
    return (alignment.hit_def, rec.query, (100.*hsp.identities/hsp.align_length), query_cov,
            rec.query_length, alignment.length, hsp.identities, hsp.align_length,
            hsp.query_start,hsp.query_end,hsp.sbjct_start,hsp.sbjct_end,
            hsp.expect,hsp.query_strand,hsp.sbjct_strand)


def iter_blast_records_to_tuples(records,
                                 hsp_filter=hsp_filter_default,
                                 hsp_filter_args={},
                                 id_extractor=r"^(\S+)"):
    """Iterate through BLAST records and yield tuples that match the criteria."""
    ##Note that this will be very different if we switch to feeding records from SeqIO module.
    ##In particular, coordinates will be Pythonic zero-based half-open and already flipped for
    ##negative strand (see Bio.SeqIO manual)
    import re
    yield hsp_filter()
    id_extractor_re = re.compile(id_extractor)
    id_extractor_f = lambda x: re.match(id_extractor_re,x).group(0)
    for rec in records:
        ## one rec per query
        ## query ID from defline
        ##rec.query
        ## length of query sequence
        ##rec.query_length
        rec.query = id_extractor_f(rec.query)
        for alignment in rec.alignments:
            ## one alignment per subject
            ## subject ID from defline
            ##alignment.hit_def
            ## length of subject sequence
            ##alignment.length
            alignment.hit_def = id_extractor_f(alignment.hit_def)
            for hsp in alignment.hsps:
                ## the next three are numbers
                ##hsp.identities
                ##hsp.positives
                ##hsp.gaps
                ## length of alignment
                ##hsp.align_length
                ## two fields that reflect strand: hsp.frame and hsp.strand. The latter
                ## is created by the parser. For some reason, it is always (None,None) while
                ## the hsp.frame shows correctly (1,-1) when hsp.sbjct_end > hsp.sbjct_start
                if rec.application.upper().endswith("BLASTN"):
                    ## We pull corresponding values from frame if strand is None
                    strand = tuple([ (x[0] if x[0] is not None else x[1]) \
                                     for x in zip(hsp.strand,
                                                  (-1 if _<0 else 1 if _>0 else 0 for _ in hsp.frame)) ])
                    hsp.strand = strand
                hsp.query_strand = hsp.strand[0]
                hsp.sbjct_strand = hsp.strand[1]
                ## and always put coords in order, similar to Bio.SeqIO. NCBI BLAST+ BLASTN reverses only subject
                _order_coords = lambda x,y: (x,y) if x <= y else (y,x)
                hsp.query_start, hsp.query_end = _order_coords(hsp.query_start, hsp.query_end)
                hsp.sbjct_start, hsp.sbjct_end = _order_coords(hsp.sbjct_start, hsp.sbjct_end)

                res = hsp_filter(rec=rec,alignment=alignment,hsp=hsp,**hsp_filter_args)
                if res:
                    yield res

def iter_blast_records_to_table(records,
                            hsp_filter=hsp_filter_default,
                            hsp_filter_args=dict(),
                            out_blast_csv="blast.txt",
                            return_pandas=False):
    """Iterate through BLAST records and create a delimited file-table for those records that match the criteria.
    Optionally, return Pandas DataFrame."""
    with open(out_blast_csv,"w") as out_bl:
        bl_writer = csv.writer(out_bl, dialect='excel-tab', lineterminator='\n')

        for rec in iter_blast_records_to_tuples(records,
            hsp_filter=hsp_filter,
            hsp_filter_args=hsp_filter_args):
            bl_writer.writerow(rec)
    if return_pandas:
        import pandas
        return pandas.read_csv(out_blast_csv,dialect="excel-tab")

def _main():
    import argh
    argh.dispatch_commands([
        makeblastdb
        ])

if __name__ == "__main__":
    _main()
