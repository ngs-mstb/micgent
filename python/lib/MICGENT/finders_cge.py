"""Implementation of a few published algorithms from the Center for Genetic Epidemiology"""
from __future__ import absolute_import
from __future__ import division
from __future__ import unicode_literals
from __future__ import print_function

from future import standard_library
standard_library.install_aliases()
from builtins import *
from past.utils import old_div
from . import util
from . import blast
from . import config
from . import annotations

from MGT.Logging import *


def hsp_filter_cge(rec=None, alignment=None, hsp=None, db_cov_min=0.6, pct_ident_min=80.):
    """Filter and convert HSP records as it is done by XXXFinder tools from CGE.
    CGE stands for Center for Genetic Epidemiology
    Queries are assumed to be the known reference sequences, and subjects -
    the unknown sequences to be classified. We swap the names of the output
    fields in this function so that unknown sequences are called queries.
    db_cov_min parameter sets cutoff on the reference sequences
    coverage.
    """
    ## PlasmidFinder paper (PMC4068535) says:
    # ```
    # Upon sequence submission, a percent identity (% ID) threshold (the percentage of
    # nucleotides that are identical between the best-matching replicon sequence in the
    # database and the corresponding sequence in the assembled sequencing data) of 100%,
    # 95%, 90%, 85%, 80%, or on down to 50% can be selected. However, it should be noted
    # that in the current form, PlasmidFinder is designed to identify replicons with at
    # least 80% nucleotide identity with those currently included in the database (see
    # Table S1 in the supplemental material) and will not adequately cover plasmid
    # diversity outside this scope.
    # ```
    # ```
    # Details on how the best match is selected are described in reference
    # PMC3318499. For a hit to be reported, it has to cover at least 60% of
    # the length of the replicon sequence in the database. Output data include
    # information on what DNA fragment (contig) was found and the position of the
    # hit within this contig. Also, information regarding the % ID, the length of
    # the hit, and the length of the replicon sequence is included in the output.
    # ```
    # The PlasmidFinder paper sends to an earlier paper on MLST (PMC3318499)
    # for the match selection criteria, which says:
    # ```
    # The best-matching MLST allele is found by calculating the length
    # score (LS) as QL - HL + G, where QL is the length of the MLST
    # allele, HL is the length of the HSP, and G is the number of gaps
    # in the HSP. The allele with the lowest LS and, secondly, with
    # the highest percentage of identity (ID) is selected as the
    # best-matching MLST allele.
    # ```
    if rec is None:
        return ("query_id", "sbjct_id", "length_score", "pct_ident", "sbjct_cov",
                "sbjct_length", "query_length", "identities", "align_length",
                'sbjct_start', 'sbjct_end', 'query_start', 'query_end', 'query_strand', 'sbjct_strand')
    nongap_align_length = hsp.align_length - hsp.gaps
    length_score = rec.query_length - nongap_align_length
    query_cov = (old_div(float(nongap_align_length), rec.query_length)) if rec.query_length > 0 else 0
    pct_ident = (100. * hsp.identities / hsp.align_length)
    if query_cov >= db_cov_min and pct_ident >= pct_ident_min:
        return (alignment.hit_def, rec.query, length_score, pct_ident, query_cov,
                rec.query_length, alignment.length, hsp.identities, hsp.align_length,
                hsp.query_start, hsp.query_end, hsp.sbjct_start, hsp.sbjct_end,
                hsp.sbjct_strand, hsp.query_strand)


def select_top_overlapping_hits_cge(bl_df,closed_end=True):
    """Use annotations.select_top_overlapping_intervals() to select CGE hits"""
    ##groupby preserves the order of rows within each group
    ##as_index=False means rows are grouped in memory, not indirectly by index
    bl_sel_df = bl_df.sort_values(by=["query_id", "length_score", "pct_ident", "sbjct_id"],
                                  ascending=[True, True, False, True]) \
        .groupby(['query_id','locus'], as_index=False) \
        .apply(lambda x: x.iloc[
                         list(annotations.select_top_overlapping_intervals(
                             x.query_start, x.query_end, closed_end=closed_end
                         )),
                         :
                         ])
    bl_sel_df.reset_index(drop=True,inplace=True)
    return bl_sel_df


def annotate_cge(seq=None,
                 seq_list=None,
                 seq_db=None,
                 seq_db_list=None,
                 db_format="fasta",
                 meta_db=None,
                 meta_db_re=None,
                 out_blast_csv="blast.txt",
                 out_annot_csv="annot.txt",
                 db_cov_min=0.6,
                 pct_ident_min=80.,
                 blast_app="blastn",
                 no_blast_cleanup=False,
                 short_reference=False,
                 db_strand="both",
                 threads=1):
    """Annotator that encapsulates common approach of FinderXXX tools from the Center for Genetic Epidemiology.
    This returns best BLASTN hits against the reference DB filtered by reference coverage and percent identity.
    Coordinates are one-based closed ranges as returned by BLAST.
    Note that more than one record per query sequence can be returned if there are non-overlapping matches
    on the query. You need to post-filter the records if the expectation is that there can be only one
    annotation per query.
    @param meta_db Delimited text file with fields sbjct_id, locus, cluster
    @param meta_db_re Python regular expression with named capture groups to extract subject metadata
    from the sbjct_id, as an alternative to providing meta_db file. Example: '^(?P<locus>[^_]+)_(?P<cluster>.+)'
    will extract locus=kdpE cluster=kdp_SCCMec from ID=kdpE_kdp_SCCMec
    """
    import pandas as pd
    blast_par = {"outfmt": 5,
                 "max_target_seqs": 1000000000,
                 "num_threads": threads}
    if db_strand != "both":
        if blast_app in ("blastn","blastx","tblastx"):
            blast_par["strand"] = db_strand
    blast_res = blast.blast(seq_query=seq_db,
                            seq_query_list=seq_db_list,
                            seq_db=seq,
                            seq_db_list=seq_list,
                            db_format=db_format,
                            app=blast_app,
                            return_records=True,
                            blast_par=blast_par,
                            short_reference=short_reference,
                            cleanup=not no_blast_cleanup)
    bl_df = blast.iter_blast_records_to_table(blast_res["records"],
                                              hsp_filter=hsp_filter_cge,
                                              hsp_filter_args=dict(db_cov_min=db_cov_min,
                                                                   pct_ident_min=pct_ident_min),
                                              return_pandas=True,
                                              out_blast_csv=out_blast_csv)
    if meta_db:
        meta_db_df = pd.read_table(meta_db,sep=None) #auto-guess sep
        sbjct_id_before = bl_db.sbjct_id.copy()
        bl_df = pd.merge(bl_df,meta_db_df,how="inner",on="sbjct_id")
        bl_df.reset_index(drop=True, inplace=True)
        assert (sbjct_id_before == bl_df.sbjct_id).all(),"Some DB records did not have matching metadata records"
    elif meta_db_re:
        bl_df = pd.concat([bl_df, bl_df.sbjct_id.str.extract(meta_db_re,expand=True)],axis=1)
        bl_df.reset_index(drop=True, inplace=True)
        assert bl_df.locus.notnull().all(), "Extracting locus field from sbjct_id failed for some fields"
    else:
        ## by default, assume that all DB sequences are alleles of a single locus
        bl_df["locus"] = "Locus1"
    assert "locus" in list(bl_df), \
        "DB metadata should at least have a 'locus' field"
    if "cluster" not in list(bl_df):
        ## bu default, each locus is its own cluster
        bl_df["cluster"] = bl_df.locus
    bl_sel_df = select_top_overlapping_hits_cge(bl_df, closed_end=True)
    if out_annot_csv:
        bl_sel_df.to_csv(out_annot_csv, index=False, sep="\t")
    return bl_sel_df


def annotate_plasmids(seq, seq_db=None,
                      db_format="fasta",
                      out_blast_csv="blast.txt",
                      out_annot_csv="annot.txt",
                      db_cov_min=0.6,
                      pct_ident_min=95.):
    """Implements the algorithm of PlasmidFinder (PMC4068535).
    Default pct_ident_min=95. is taken from the setting on the
    Web server (2016-04-29).
    Coordinates are one-based closed ranges as returned by BLAST"""
    if not seq_db:
        seq_db = config.load_config(pkg="MICGENT")["data"]["cge"]["plasmid_db"]["seq"]
    return annotate_cge(seq=seq, seq_db=seq_db,
                        db_format=db_format,
                        out_blast_csv=out_blast_csv,
                        out_annot_csv=out_annot_csv,
                        db_cov_min=db_cov_min,
                        pct_ident_min=pct_ident_min)


def annotate_resistance_genes(seq, seq_db=None,
                              db_format="fasta",
                              out_blast_csv="blast.txt",
                              out_annot_csv="annot.txt",
                              db_cov_min=0.4,
                              pct_ident_min=100.):
    """Implements the algorithm of ResFinder (PMC3468078)
    Coordinates are one-based closed ranges as returned by BLAST"""
    # All genes from the ResFinder database were BLASTed against the assembled genome, and
    # the best-matching genes were given as output. For a gene to be reported, it has to
    # cover at least 2/5 of the length of the resistance gene in the database. The best-matching
    # genes were identified as previously. It is possible to select a % identity (ID) threshold
    # (the percentage of nucleotides that are identical between the best-matching resistance gene
    # in the database and the corresponding sequence in the genome). The default ID is 100%.
    if not seq_db:
        seq_db = config.load_config(pkg="MICGENT")["data"]["cge"]["res_db"]["seq"]

    bl_df = annotate_cge(seq=seq, seq_db=seq_db,
                         db_format=db_format,
                         out_blast_csv=out_blast_csv,
                         out_annot_csv=None,
                         db_cov_min=db_cov_min,
                         pct_ident_min=pct_ident_min)
    bl_df["res_class"] = bl_df.sbjct_id.str.rsplit("_", 1).str[0]
    if out_annot_csv:
        bl_df.to_csv(out_annot_csv, index=False, sep="\t")
    return bl_df


## import package module and add argh entry points

def _main():
    import argh
    argh.dispatch_commands([
        annotate_plasmids,
        annotate_resistance_genes,
        annotate_cge
    ])


if __name__ == "__main__":
    _main()
