"""Format converters"""
from __future__ import absolute_import
from __future__ import unicode_literals
from __future__ import print_function
from __future__ import division

from future import standard_library
standard_library.install_aliases()
from builtins import range
from builtins import *
from . import util
from MGT.Logging import *

def _gto_feature_coords(loc):
    """Return coordinates from GTO location field as half-open zero based range"""
    #GTO native converter does this: [804,-,237] -> [568,804,-]
    #when converting to GFF (which is supposed to be 1-based inclusive range).
    #We return [567,804,-] as in Python or BED 0-based semi-open range.
    start = int(loc[1])
    length = int(loc[3])
    strand = loc[3].strip()
    if strand == "-":
        end = start
        start = end - length
    else:
        start -= 1
        end = start + length
    return (loc[0],start,end,strand)


def gto_to_gff(gto_file,gff_file):
    """Convert RASTtk GTO JSON file into GFF3 file.
    Compared to a converter provided by RASTtk,
    this one includes more information"""
    import MGT.GFF as _gff
    import json
    with open(inp_file,"r") as inp:
        js = json.load(inp)
    #js["contigs"] is a dictionary of contig sequences
    #
    with open(gff_file,"w") as gff_out:
        _gff.GFF3Header().write(gff_out)
        for x in js["contigs"]:
            _gff.GFF3SeqRegion(x["id"],1,len(x["dna"])).write(gff_out)
        for x in js["features"]:
            loc = x["location"]
            assert len(loc) == 1, "I do not know how to treat multi-location feature: {}"
            loc = _gto_feature_coords(loc[1])
            attribs=dict(
                ID=x["id"],
                Name=x["function"]
            )
            _gff.GFF3Record(seqid=loc[0],
                            source="RAST2",
                            type=x["type"],
                            start=loc[1]+1,
                            end=loc[2],
                            strand=loc[3],
                            attribs=attribs).write(gff_out)
        _gff.GFF3Fasta().write(gff_out)
        for x in js["contigs"]:
            _gff.GFF3Fasta().write(gff_out)
            import Bio
            rec = Bio.SeqRecord.SeqRecord(Bio.Seq.Seq(x["dna"], Bio.Alphabet.generic_dna),
                 id=x["id"],
                 description="")
            SeqIO.write([rec], gff_out, "fasta")


def genbank_to_gff(gff_file,fasta_file=None,genbank_files=None,
                   genbank_files_list=None,include_seq=False,
                   rec_id_field="id"):
    from Bio import SeqIO
    from BCBio import GFF
    out_fasta = None
    try:
        if fasta_file:
            out_fasta = open(fasta_file,"w")
            include_seq = True
        with open(gff_file, "w") as out_handle:
            for genbank_file in util.glob_files(files_globs=genbank_files,
                                         files_globs_list=genbank_files_list):
                log.info("Processing sequences from {}".format(genbank_file))
                GFF.write(SeqIO.parse(genbank_file, "genbank"), out_handle,
                          include_fasta=include_seq,
                          out_handle_fasta=out_fasta,
                          rec_id_field=rec_id_field)
    finally:
        if out_fasta:
            out_fasta.close()


def label_gff_features_with_annotations(gff_file,
                                        annot_file,
                                        feature_type="CDS",
                                        annot_join="interval",
                                        annot_sep="\t"):
    """Join features in GFF3 file with annotations.
    Add annotations as attributes of GFF features.
    The way to join with the annotation file must be provided in annot_join, which
    can be either keyed by nucleotide intervals in the first three fields like GFF itself, or
    by feature IDs (if annotation was done for protein sequences - not
    yet implemented)
    """
    import petl
    import petlx
    feat = petlx.bio.gff3.fromgff3(gff_file)
    ann = petl.io.csv.fromcsv(annot_file,delimiter=annot_sep)
    ann = petl.convert(ann, ('query_start', 'query_end','query_strand'), int)
    jn = petl.transform.intervals.intervaljoin(feat, ann, lstart="start", lstop="end", rstart="query_start",
                                          rstop="query_end", lkey="seqid", rkey="query_id", include_stop=True)


def petl_fix_gff_coord(feat):
    """For petl tables loaded from GFF3 files, make the coords and strand standard-compliant.
    For example, as of 2016-05-02, GFF3 files produces by RASTtk might have start > end with positive strand
    (for phage finder annotations, and additionally ID that is not unique within the GFF3 file).
    """
    import petl

    ## result either `-` or `.`
    _lk = {'+': '-', '-': '-', ".": "."}
    feat = petl.convert(feat,
                      {
                          'start': lambda v, row: row.end,
                          'end': lambda v, row: row.start,
                          'strand' : lambda v, row, _lk=_lk: lk[v]
                      },
                      where = lambda row: row.start > row.end,
                      pass_row = True
                      )
    return feat

import contextlib
@contextlib.contextmanager
def petl_opened_file_source(filename,mode,**kw):
    class _petl_opened_file_source_class(file):

        def open(self,*l,**kw):
            return self

        def close(self):
            pass

        def really_close(self):
            super(_petl_opened_file_source_class,self).close()

    f = None
    try:
        f = _petl_opened_file_source_class(filename,mode=mode,**kw)
        yield f
    finally:
        if f:
            f.really_close()

import json

def encode_json_fasta_defline(descr):
    return json.dumps(descr, ensure_ascii=True, indent=None)

def decode_json_fasta_defline(descr):
    return json.loads(descr)

def join_annotations_with_gff_features(annot_file,
                                       annot_file_out,
                                       feature_type="CDS",
                                       annot_join="interval",
                                       annot_sep="\t",
                                       gff_files=None,
                                       gff_files_list=None,
                                       max_overlap_only=True):
    """Join features in GFF3 file with annotations.
    Add annotations as attributes of GFF features.
    The way to join with the annotation file must be provided in annot_join, which
    can be either keyed by nucleotide intervals in the first three fields like GFF itself, or
    by feature IDs (if annotation was done for protein sequences - not
    yet implemented).
    Coordinates and strand of each annotation will be replaced with those of overlapping feature,
    if any (can result in more than one record per annotation).
    """
    import petl
    import petlx

    ann_all = petl.io.csv.fromcsv(annot_file, delimiter=annot_sep)
    ann_all = petl.convert(ann_all, ('query_start', 'query_end', 'query_strand'), int)
    ann_all = petl.addcolumn(ann_all,'ann_rec_ind',range(ann_all.nrows()))

    with petl_opened_file_source(annot_file_out,"w") as annot_out:
        for i_inp, gff_file in enumerate(util.glob_files(files_globs=gff_files,
                                                         files_globs_list=gff_files_list)):
            log.info("Working on feature file {}".format(gff_file))
            feat = petlx.bio.gff3.fromgff3(gff_file)
            feat = petl_fix_gff_coord(feat)

            feat_seqid_set = set(feat["seqid"])

            if feature_type:
                feat = petl.selecteq(feat,'type',feature_type)

            ann = petl.selectin(ann_all, 'query_id', feat_seqid_set)

            ## somehow we get many ORFs in GFFs (and Genbank files) from both RASTtk and ClovR where one
            ## ORFs ends at the start position of another ORF (and the BLAST match starts at the start of the
            ## second ORF).
            jn = petl.transform.intervals.intervalleftjoin(ann, feat, rstart="start", rstop="end", lstart="query_start",
                                                  lstop="query_end", rkey="seqid", lkey="query_id",
                                                           rprefix="feat_")
            jn = petl.addfield(jn,"overlap_len",
                               lambda rec: (min(rec['end'],rec['query_end']) - max(rec['start'],rec['query_start']) + 1) \
                                       if rec['start'] is not None else 0)
            if max_overlap_only:
                jn = petl.groupselectmax(jn,key="ann_rec_ind",value="overlap_len")
            _strand_conv = {'+':1, '-':-1, '.':0}
            jn = petl.convert(jn,
                              {
                                  'query_start' : lambda v,row: row.start if row.start is not None else v,
                                  'query_end': lambda v,row: row.end if row.end is not None else v,
                                  'query_strand': lambda v,row,_strand_conv=_strand_conv: _strand_conv[row.strand] \
                                      if row.strand is not None else row.query_strand
                              },
                              pass_row=True
                              )
            if i_inp == 0:
                out_func = petl.io.csv.tocsv
            else:
                out_func = petl.io.csv.appendcsv
            out_func(jn, annot_out, delimiter=annot_sep)

def annotations_to_gff(annot_file,
                       gff_file,
                       feature_type="CDS",
                       feature_source="MICGENT",
                       name_field="sbjct_id",
                       annot_sep="\t"):
    """Create GFF3 file out of delimited annotation file.
    This assumes one-based closed ranges in annotation file"""
    import MGT.GFF as _gff
    import MGT.UUID as _uuid
    import pandas as pd
    df = pd.read_table(annot_file,sep=annot_sep)
    df.sort_values(by=["query_id","query_start","query_end","sbjct_id","sbjct_start"],inplace=True)
    cols = list(df)
    cols_attrib = [ _ for _ in cols if _ not in ("query_id","query_start","query_end","query_strand") ]
    with open(gff_file,"w") as gff_out:
        _gff.GFF3Header().write(gff_out)
        for g in df.groupby("query_id"):
            query_id = g[0]
            recs = g[1]

            _gff.GFF3SeqRegion(query_id,1,recs.iloc[0].query_length).write(gff_out)

            ##iterrows converts each row to a series and destroys types
            for i in range(recs.shape[0]):

                rec = recs.iloc[i]

                attribs=dict(
                    ID=_uuid.genId(),
                    Name=rec[name_field]
                )
                attribs.update({name:rec[name] for name in cols_attrib})
                # this class is zero based
                _gff.GFF3Record(seqid=query_id,
                                source=feature_source,
                                type=feature_type,
                                start=rec.query_start-1,
                                end=rec.query_end,
                                strand=rec.query_strand,
                                attribs=attribs).write(gff_out)


def extract_padded_range(seq_rec,start,end,pad_left,pad_right,rev_comp,seq_rec_id=None,
                         annot=None,trim_annot=True):
    """Extract a range from a sequence record padding the range and reverse complement if required.
    The input range is a closed, one-based range, to be consistent with the annotation.
    @return dict(seq=seq rec out, start, end=range relative to seq rec out), coords are one-based"""
    #assert pad_left >= 0 and pad_right >= 0
    start -= 1
    seq_len = len(seq_rec)
    start_do = min(max(start - pad_left,0),seq_len)
    end_do = max(min(end + pad_right,seq_len),0)
    assert start_do <= end_do
    seq_rec_out = seq_rec[start_do:end_do]
    #import ipdb; ipdb.set_trace()
    start_out = start - start_do
    end_out = end - start_do
    seq_len_out = len(seq_rec_out)
    if trim_annot:
        start_out = min(max(start_out,0),seq_len_out)
        end_out = min(max(end_out, 0),seq_len_out)
    if not seq_rec_id:
        seq_rec_id = seq_rec.id
    if annot is not None:
        assert (annot.query_id==seq_rec.id).all(),"All annotations must be defined on the target sequence"
        annot = annot.copy()
        annot.query_start -= 1
        annot.query_start -= start_do
        #assert (annot.query_start >= 0).all() # or we can use .clip()
        annot.query_end -= start_do
        annot.query_id = seq_rec_id
        annot.query_length = seq_len_out
    if rev_comp:
        seq_rec_out = seq_rec_out.reverse_complement()
        start_out_pre = start_out
        start_out = seq_len_out - end_out
        end_out = seq_len_out - start_out_pre
        if annot is not None:
            start_out_pre = annot.query_start.copy()
            annot.query_start = seq_len_out - annot.query_end
            annot.query_end = seq_len_out - start_out_pre
            annot.query_strand *= -1
    seq_rec_out.id = seq_rec_id
    seq_rec_out.description = seq_rec.description
    if annot is not None:
        annot.query_start += 1
    return (seq_rec_out,
            start_out+1,
            end_out,
            annot)

def extract_sequence_for_annotations(annot_file,
                                     fasta_file_out,
                                     pad=0,
                                     pad_left=None,
                                     pad_right=None,
                                     pad_left_origin="left",
                                     pad_right_origin="right",
                                     annot_file_out=None,
                                     split_by=None,
                                     stats_file=None,
                                     annot_sep='\t',
                                     fasta_files=None,
                                     fasta_files_list=None,
                                     omit_cluster_id=False,
                                     keep_seq_non_matching=False):
    """Extract nucleotide sequences corresponding to annotations.
    This assumes one-based closed ranges in annotation file.
    This will reverse-complement the sequence if the annotation is given on
    the negative strand."""
    from Bio import SeqIO
    import pandas as pd
    import re
    import json
    import csv

    if pad_left is None:
        pad_left = pad
    else:
        pad_left = int(pad_left)
    if pad_right is None:
        pad_right = pad
    else:
        pad_right = int(pad_right)

    if annot_file_out is None:
        annot_file_out = util.rm_fname_extension(fasta_file_out)+".txt"

    ann_all = pd.read_table(annot_file,sep=annot_sep)
    ann_fac = ann_all.groupby("query_id")
    ann_by = {}
    i_rec_out = 0
    ann_writer = None
    with open(fasta_file_out,"w") as fasta_out,\
        open(annot_file_out,"w") as annot_out:
        for fname in util.glob_files(files_globs=fasta_files,
                                     files_globs_list=fasta_files_list):
            log.info("Reading sequences from {}".format(fname))
            for rec_cont in SeqIO.parse(fname, "fasta"):
                if rec_cont.id in ann_fac.groups:
                    ann_cont = ann_fac.get_group(rec_cont.id)
                    ann_clust = ann_cont.groupby("cluster")
                    if omit_cluster_id:
                        assert len(ann_clust) == 1
                    for cluster,rows in ann_clust:
                        start = rows.query_start.min()
                        end = rows.query_end.max()
                        start_extract=start.copy()
                        end_extract=end.copy()
                        if pad_left_origin=="right":
                            start_extract=end
                        if pad_right_origin=="left":
                            end_extract=start
                        rev_comp = (rows.query_strand.sum()<0)
                        assert len(cluster.split())==1, "Spaces are not allowed in cluster ID"
                        rec_ann,start_ann,end_ann,rows_ann = extract_padded_range(rec_cont,
                                                                 start=start_extract,
                                                                 end=end_extract,
                                                                 pad_left=pad_left,
                                                                 pad_right=pad_right,
                                                                 rev_comp=rev_comp,
                                                                 seq_rec_id="{}_{}".format(rec_cont.id,cluster) \
                                                                                  if not omit_cluster_id \
                                                                     else rec_cont.id,
                                                                 annot=rows)
                        #rec_ann.description = encode_json_fasta_defline([row_ann,row])
                        rows_ann.sort_values(by=["query_id", "query_start", "query_end", "sbjct_id", "sbjct_start"],
                                       inplace=True)
                        SeqIO.write([rec_ann],fasta_out,"fasta")
                        if i_rec_out==0:
                            ann_writer = csv.writer(annot_out,
                                                        dialect='excel-tab',
                                                        lineterminator='\n',
                                                        delimiter=annot_sep)
                            ann_writer.writerow(list(rows_ann))
                        ann_writer.writerows(rows_ann.itertuples(index=False))
                        if split_by:
                            split_by_val = rows_ann[split_by]
                            if not (split_by_val==split_by_val.iloc[0]).all():
                                raise ValueError("Split attribute has more than one value within cluster {}".format(cluster))
                            ann_by.setdefault(split_by_val.iloc[0],[]).append(rec_ann.id)
                        i_rec_out += 1
                elif keep_seq_non_matching:
                    SeqIO.write([rec_cont], fasta_out, "fasta")

    if split_by:
        assert stats_file, "If splitting, you need to provide a name for the output stats file"
        stats = dict(splits=dict(key=split_by,files=dict()))
        files = stats["splits"]["files"]
        seq_dict = SeqIO.index(fasta_file_out, "fasta")
        for key, rec_anns in list(ann_by.items()):
            key_file = re.sub(r"\s+","-",key)
            key_file = re.sub(r"/+","_",key_file)
            fasta_split_file = "{}_{}.fasta".format(fasta_file_out,key_file)
            with open(fasta_split_file,"w") as fasta_out:
                for rec_id in rec_anns:
                    SeqIO.write([seq_dict[rec_id]],fasta_out,"fasta")
            files[key] = dict(seq_file=fasta_split_file)
        with open(stats_file,"w") as out:
            json.dump(stats,out,indent=4)

def cat_lines(files=None,
          files_list=None,
              remove_all_empty=False,
              keep_last_empty=False,
              skip_header=0,
              keep_first_header=True):
    """Concatenate lines from multiple files.
    This is similar to Unix `cat -s` command, but works
    consistently across platforms, removes all empty lines
    between files and takes globs without expanding
    them on the command line in order not to be limited
    by the OS command line length limit.

    Primary use case is concatenating FASTA files, any of
    which might or might not have eol at the end of last
    record.
    @param remove_all_empty remove also empty lines in the
    middle of each file
    @param eol symbol at the end of each file becomes extra
    eol
    @result writes concatenated files to stdout
    """
    import sys
    files_inp = list(util.glob_files(files_globs=files,
                                     files_globs_list=files_list))
    out = sys.stdout
    last_line = ""
    for i_file,f_inp in enumerate(files_inp):
        with open(f_inp,"r") as inp:
            for i_line,line in enumerate(inp):
                if i_line >= skip_header or (keep_first_header and i_file==0):
                    out.write(last_line)
                    last_line = line
                    if remove_all_empty:
                        last_line = last_line.rstrip("\n")
                        if last_line or keep_last_empty:
                            last_line += "\n"
            else:
                last_line = last_line.rstrip("\n")
                if last_line or keep_last_empty:
                    last_line += "\n"

    out.write(last_line)

## import package module and add argh entry points

def _main():
    import argh
    argh.dispatch_commands([
        join_annotations_with_gff_features,
        annotations_to_gff,
        extract_sequence_for_annotations,
        genbank_to_gff,
        cat_lines
        ])

if __name__ == "__main__":
    _main()
