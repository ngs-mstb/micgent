# Copyright Medimmune LLC, 2016
# Copyright J. Craig Venter Institute, 2013-2015
#
# The creation of this software was supported by the DARPA Prophecy
# program.
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

"""Methods for working with alignments and sequences"""
from __future__ import print_function
from __future__ import absolute_import
from __future__ import division
from __future__ import unicode_literals

from future import standard_library
standard_library.install_aliases()
from builtins import zip
from builtins import str
from builtins import range
from builtins import *
from builtins import object
from past.utils import old_div
from .bio_imports import *
from .std_imports import *
from . import util

import numpy as np
import pandas as pd
import collections


def get_dna_ambiguous_map(to_amb=True):
    from Bio.Data import IUPACData
    ##"B":"CGT"
    to = IUPACData.ambiguous_dna_values
    if to_amb:
        return to
    else:
        m = {}
        for a,na in list(to.items()):
            na = "".join(sorted(na))
            ##both X and N map to CGTA, use only N in reversed map
            if a == "X":
                a = "N"
            m[na] = a
        return m

def str_as_np(s):
    ##dtype=np.character will make something like '|S347', not
    ##'|S1' that we need.
    return np.array(s,'c')

def ali_as_np(ali,order="C"):
    """Use order='F' when mostly slicing along columns"""
    ali_np = np.array([str(rec.seq) for rec in ali], 'c', order=order, ndmin=2)
    if len(ali_np.shape)< 2:
        ## looks like a bug caused by Numpy assignment operator overloading:
        ## if I simply have r = r[:,None], then r still has shape==(1,).
        ## coding it as a work-around
        ali_np = ali_np[:,None]
    assert ali_np.shape == (len(ali),ali.get_alignment_length())
    return ali_np

def shift_2d_col(x,n,fill=0):
    x = np.roll(x,n,axis=1)
    if n > 0:
        x[:,:n] = fill
    else:
        x[:,n:] = fill
    return x

def shift_1d(x,n,fill=0):
    x = np.roll(x,n)
    if n > 0:
        x[:n] = fill
    else:
        x[n:] = fill
    return x

def pad_index_1d(x,pad,min_val,max_val):
    return np.unique(np.sort(np.clip(np.concatenate([x,x-pad,x+pad]),min_val,max_val)))

def pos_mask_from_ind(pos_ind,len_seq):
    if pos_ind is None:
        pos_mask = np.ones(len_seq,dtype=bool)
    else:
        pos_mask = np.zeros(len_seq,dtype=bool)
        pos_ind = np.asarray(pos_ind,dtype=int)
        pos_mask[pos_ind] = True
    return pos_mask

def ali_coord_map(ali_np,gap_symb="-"):
    """Create coordinate mappings between alignment and its sequences.
    @param ali_np - 2d numpy characater array of alignment
    @return dict with mappings for alignment coords to sequence coords (2d array)
    called ali_to_seq
    and sequence coords to alignment coords (list of 1d arrays)
    called seq_to_ali. ali_to_seq will have non-increasing numbers at
    positions corresponding to gaps. If the alignment starts with a gap,
    the corresponding row in ali_to_seq will start with a zero, the same as if
    the alignment starts with a non-gap. This ensures that all coordinates are
    always valid positions in the sequence, to be used in operations like "lift-over".
    To remove gap positions, you will need to use something like
    ali_to_seq[np.where(ali_np != gap_symb)].
    All coords are zero-based.
    If you want to get original sequence i, use ali_np[i,self["seq_to_ali"][i]].
    To make alignment string out of original sequence, use ali_np[i,self["seq_to_ali"][i]] = s
    """
    #start of sequence is first >= 0 value:
    ali_to_seq = np.arange(ali_np.shape[1]) - np.cumsum(ali_np == gap_symb,axis=1)
    #get 1s at each non-gap position - pad with -1 at left and find diff with previous element
    non_gaps = np.diff(np.lib.pad(ali_to_seq,((0,0),(1,0)),"constant",constant_values=-1),axis=1)
    seq_to_ali = [ np.nonzero(non_gaps[i,:])[0] for i in range(non_gaps.shape[0]) ]
    ali_to_seq[np.where(ali_to_seq < 0)] = 0
    return dict(ali_to_seq=ali_to_seq,seq_to_ali=seq_to_ali)

def seq_id_to_ind(seqs):
    """Iterate through Bio.alignment or some other sequence set and return dict mapping IDs to index"""
    m = {}
    for i,a in enumerate(seqs):
        assert a.id not in m, "Duplicate sequence ID: {}".format(a.id)
        m[a.id] = i
    return m

def _ind_to_np(ind):
    if not util.is_sequence(ind):
        ind = [ ind ]
    return np.asarray(ind)

def _ind_to_np2(ind1,ind2):
    return (_ind_to_np(ind1),
            _ind_to_np(ind2))

def ali_np_ident_matrix(ali_np,ind1,ind2):
    """Return len(ind1) x len(ind2) matrix of identity counts between two partitions of sequences in the alignment"""
    ind1, ind2 = _ind_to_np2(ind1,ind2)
    res = np.zeros((len(ind1),len(ind2)))
    x2 = ali_np[ind2]
    for (i,ind) in enumerate(ind1):
        res[i,:] = (ali_np[ind] == x2).sum(1)
    return res

## pairwise alignment length is: extract alignment for the pair;
## cut pair alignment to overlap of start, end ranges
## and compute the number of not-all-gap columns.
## That gives a matrix. Coverage is computed as two
## matrices by dividing the pairwise alignment length matrix
## by the lengths of row and column sequences respectively

def ali_np_pairwise_length(ali_np,ali_range,ind1,ind2,gap_symb="-"):
    ind1, ind2 = _ind_to_np2(ind1,ind2)
    ali_len = np.zeros((len(ind1), len(ind2)))
    ##TODO: express outer loop in array operations for speed
    for (i1, in1) in enumerate(ind1):
        for (i2, in2) in enumerate(ind2):
            a12 = ali_np[(in1,in2),
            min(ali_range[in1,0],ali_range[in2,0]):
            max(ali_range[in1, 1], ali_range[in2, 1])]
            a12_len = ((a12==gap_symb).sum(0) == 2)
            ali_len[i1,i2] = a12_len
    return ali_len

def ali_np_pairwise_coverage(ali_len,seq_len,ind1,ind2):
    ali_range1 = ali_range[ind1]
    ali_range2 = ali_range[ind2]
    ali_len1 = ali_range1[:, 1] - ali_range1[:, 0]
    ali_len2 = ali_range2[:, 1] - ali_range2[:, 0]
    cov1 = old_div(ali_len, ali_)
    return res


class ali_mapper(object):
    """Class that opens an MSA file and also caches some derived data
    structures"""

    def __init__(self,
            msa_file,
            msa_format="fasta",
            ali_np_order="F",
            gap_symb="-"):
        if util.is_string(msa_file):
            self.ali = AlignIO.read(msa_file,msa_format)
        else:
            ## must be already an MSA object
            self.ali = msa_file
        self.ali_np = ali_as_np(self.ali,order=ali_np_order)
        self.coord_map = ali_coord_map(self.ali_np,gap_symb=gap_symb)
        self.id_to_ind = seq_id_to_ind(self.ali)
        self.gap_symb = gap_symb

    def get_bio_ali(self):
        return self.ali

    def dim(self):
        return self.ali_np.shape

    def len_seq_by_ind(self,ind=None):
        seq_to_ali = self.coord_map["seq_to_ali"]
        if ind is None:
            ind = np.arange(len(seq_to_ali))
        if not util.is_sequence(ind):
            return len(seq_to_ali[ind])
        else:
            return np.fromiter((len(seq_to_ali[_]) for _ in ind),np.int)

    def len_seq_by_id(self,id):
        return self.len_seq_by_ind(self.get_id_to_ind(id))

    def get_id_to_ind(self,id=None):
        if id is None:
            return self.id_to_ind
        id_to_ind = self.id_to_ind
        if util.is_sequence(id):
            return np.fromiter((id_to_ind[_] for _ in id),np.int)
        assert util.is_string(id)
        return id_to_ind[id]

    def get_seq_by_ind(self,ind,a=None):
        """Translate string in alignment coords into string
        coords from row ind of the alignment.
        @param ind index of alignment row
        @param a 1d array alignment string to translate. If None, get sequence
        from the alignment row 
        """
        ind_tr = self.coord_map["seq_to_ali"][ind]
        if a is None:
            return self.ali_np[ind,ind_tr]
        else:
            return a[ind_tr]

    def get_seq_by_id(self,id,a=None):
        return self.get_seq_by_ind(self.id_to_ind[id],a=a)

    def get_ali_range(self):
        """Return 2d array with start and end coordinates of each sequence in the alignment"""
        return np.asarray([(x[0], x[-1]) for x in self.coord_map["seq_to_ali"]])

    def get_seq_to_ali(self):
        """Return a list of 1D arrays N_seq x Len_seq with position of each ungapped sequence base in the alignment.
        @see ali_coord_map for the details"""
        return self.coord_map["seq_to_ali"]

    def get_ali_to_seq(self):
        """Return 2D array N_seq x N_ali_cols with position of each alignment column in each ungapped sequence.
        @see ali_coord_map for the details"""
        return self.coord_map["ali_to_seq"]

    def get_ali_np(self):
        """Return a 2D array N_seq x N_ali_cols with the alignment.
        @see ali_coord_map for the details"""
        return self.ali_np

    def get_ali_by_ind(self,ind,s=None,gap_symb="-"):
        """Take 1d array corresponding to original sequence
        from index ind in the alignment and generate an array
        in alignment coordinates. Input array can be sequence or
        something else (use proper gap_symb for gaps).
        @param s original sequence as numpy 1d array. If None,
        return sequence from alignment
        @param ind use coords for this row position in alignment
        @param gap_symb use this to fill gaps
        @return 1d array corresponding to alignment coords and dtype
        taken from s
        """
        if s is None:
            return self.ali_np[ind,:]
        r = np.ones(self.ali_np.shape[1],dtype=s.dtype)
        r[:] = gap_symb
        r[self.coord_map["seq_to_ali"][ind]] = s
        return r
    
    def get_ali_by_id(self,id,s=None,gap_symb="-"):
        return self.get_ali_by_ind(ind=self.id_to_ind[id],s=s,gap_symb=gap_symb)

    def ali_coords_by_ind(self,ind,seq_coords):
        """Convert sequence coords into alignment coords for alignment row
        ind"""
        return self.coord_map["seq_to_ali"][ind][seq_coords]

    def ali_coords_by_id(self,id,seq_coords):
        return self.ali_coords_by_ind(
                ind=self.id_to_ind[id],
                seq_coords=seq_coords
                )

    def seq_coords_by_ind(self,ind,ali_coords,remove_gaps=True):
        """Convert alignment coords into sequence coords for alignment row
        ind
        @param remove_gaps If True, the return array will have negative
        values (corresponding to gaps in alignment) removed. This might
        change the length of the output."""
        r = self.coord_map["ali_to_seq"][ind,ali_coords]
        if remove_gaps:
            r = r[np.where(ali_np[ind,ali_coords] != self.gap_symb)]
        return r

    def seq_coords_by_id(self,id,ali_coords):
        return self.seq_coords_by_ind(
                ind=self.id_to_ind[id],
                ali_coords=ali_coords
                )

    def seq_to_seq_coords_by_ind(self,ind_from,ind_to,seq_coords,remove_gaps=True):
        """Convert coords from one sequence to another sequence based on alignment.
        @param remove_gaps If True, the return array will have negative
        values (corresponding to gaps in alignment) removed. This might
        change the length of the output."""

        ali_coords = self.ali_coords_by_ind(ind_from,seq_coords)
        return self.seq_coords_by_ind(
                ind_to,
                ali_coords=ali_coords,
                remove_gaps=remove_gaps
                )

    def ident_matrix_by_ind(self,ind1,ind2):
        """Return len(ind1) x len(ind2) matrix of identity counts between two partitions of sequences in the alignment"""
        return ali_np_ident_matrix(self.ali_np,ind1,ind2)

    def ident_matrix_by_id(self,id1,id2):
        """Return len(id1) x len(id2) matrix of identity counts between two partitions of sequences in the alignment"""
        return ali_np_ident_matrix(self.ali_np,self.get_id_to_ind(id1),self.get_id_to_ind(id2))

    def seq_to_seq_coords_by_id(self,id_from,id_to,seq_coords,remove_gaps=True):
        return self.seq_to_seq_coords_by_ind(
                ind_from=self.id_to_ind[id_from],
                ind_to=self.id_to_ind[id_to],
                seq_coords=seq_coords,
                remove_gaps=remove_gaps
                )


def join_alis(alis,id_seq_keys,gap_symbol="-",error_on_loss=False):
    """Join two or more MSAs using one of the sequences in each as key.
    @param alis iter of MSAs
    @param id_seq_keys iter of IDs of key sequence in alis
    @return MSA with sum(num_rows(ali1) - 1) + 1
    (key sequence included only once)
    @pre key sequence must be identical in all alignments
    @param error_on_loss When the key sequence in alignments second and up has 
    gaps, this method will put gaps in all other sequences merged
    from those alignments (because sequences are first mapped to
    the ungapped key sequence and then to the target alignment. Set this to True
    to raise exception in such cases
    """
    for (i_ali,(ali,id_seq_key)) in enumerate(zip(alis,id_seq_keys)):
        id_seq_to_ind = seq_id_to_ind(ali)
        ali_np = ali_as_np(ali)
        coord_map = ali_coord_map(ali_np)
        key_ind = id_seq_to_ind[id_seq_key]
        ali_to_key = coord_map["seq_to_ali"][key_ind]
        if i_ali == 0:
            ali_r = ali[:] # returns new ali
            id_seq_to_ind_r = id_seq_to_ind
            ali_np_r = ali_np
            ali_to_key_r = ali_to_key
            id_seq_key_r = id_seq_key
            coord_map_r = coord_map
            key_ind_r = key_ind
        else:
            if error_on_loss:
                assert not gap_symbol in ali_np[key_ind],\
                        "Gaps in key sequence cause loss of sequence content"
            seqs_in_key = ali_np[:,ali_to_key]
            assert (ali_np[key_ind,ali_to_key] == ali_np_r[key_ind_r,ali_to_key_r]).all(),\
                    "Sequences that are used as merge keys must be identical"
            ali_in_key = np.ones((len(seqs_in_key),ali_np_r.shape[1]),"|S1")
            ali_in_key[:,:] = gap_symbol
            ali_in_key[:,ali_to_key_r] = seqs_in_key
            for i_seq, seq_np in enumerate(ali_in_key):
                if i_seq != key_ind:
                    seq_rec = ali[i_seq]
                    ali_r.append(
                        SeqRecord(
                                id=seq_rec.id,
                                description=seq_rec.description,
                                seq = Seq(seq_np.tostring(),seq_rec.seq.alphabet)
                                )
                    )
    return ali_r


def join_two_alis(ali_in1,ali_in2,id_seq_key1,id_seq_key2,ali_out,gap_symbol="-",msa_format="fasta"):
    ali1 = AlignIO.read(ali_in1,msa_format)
    ali2 = AlignIO.read(ali_in2,msa_format)
    ali = join_alis(alis=(ali1,ali2),id_seq_keys=(id_seq_key1,id_seq_key2),gap_symbol=gap_symbol)
    AlignIO.write(ali,ali_out,msa_format)

def _test_join_alis():
    seqs = [
            "---AAAAAA",
            "-CC--CCCC",
            "-LLLL--LL",
            "-TTTT--TT",
            "-GGGG--GG"
            ]
    ali = MultipleSeqAlignment(
            [
            SeqRecord(
                id="{}".format(i),
                seq=Seq(seqs[i],IUPAC.protein)
                ) \
            for i in range(len(seqs))
            ]
            )
    ali_r = join_alis(
            alis=(ali[:3],ali[2:]),
            id_seq_keys=("2","2"),
            gap_symbol="-")
    print("ali =\n", ali)
    print("ali_r =\n", ali_r)

    assert (ali_as_np(ali) == ali_as_np(ali_r)).all()

    seqs2 = [
            "LLLLLL",
            "AAAAAA",
            "CCCCCC",
            "TTTTTT",
            "GGGGGG"
            ]
    ali2 = MultipleSeqAlignment(
            [
            SeqRecord(
                id="{}".format(i),
                seq=Seq(seqs2[i],IUPAC.protein)
                ) \
            for i in range(len(seqs2))
            ]
            )
    ali_r = join_alis(
            alis=(ali,ali2),
            id_seq_keys=("2","0"),
            gap_symbol="-",
            error_on_loss=True)
    print("ali =\n", ali)
    print("ali2 =\n", ali2)
    print("ali_r =\n", ali_r)

class pdb_seq_spec(collections.namedtuple("pdb_seq_spec","chain resn resi ins")):
    __slots__ = ()

    @property
    def pos_resfile(self):
        return (str(self.resi)+str(self.ins)).rjust(4) + str(self.chain if self.chain.strip() else "_").rjust(3)
        return "{}{} {}".format(
                self.resi,
                self.ins,
                self.chain if self.chain.strip() else "_"
                )

    @property
    def resnum(self):
        assert not self.ins,"Rosetta resnum spec does not allow for insertion codes"
        ##Rosetta docs do not say how to handle empty chain ID, I am assuming here the
        ##same way as for res file ('' -> '_')
        return "{}{}".format(
                self.resi,
                self.chain if self.chain.strip() else "_"
                )
    
    @property
    def resn_one(self):
        return three_to_one(self.resn)

class pdb_seqs(dict):
    """Class that is returned from pdb_sequence()"""
    
    def get_chain(self,id):
        return self["chains_map"][id]
    
    def get_seq(self,id,format="str"):
        seq = self.get_chain(id)["seq_rec"]
        if format == "str":
            return str(seq.seq)
        elif format == "np":
            return str_as_np(str(seq.seq))
        elif format == "SeqRecord":
            return seq

def pdb_sequence(pdb_file,id=None,method="order"):
    from Bio.PDB import PDBParser, CaPPBuilder
    from Bio.PDB.Polypeptide import three_to_one
    if id is None:
        id = util.make_id_from_file_name(pdb_file)
    parser = PDBParser()
    structure = parser.get_structure(id, pdb_file)
    seq_chains = []
    for chain in structure.get_chains():
        id_chain = chain.get_id()
        if method == "distance":
            ppb=CaPPBuilder()
            seq = sum((pp.get_sequence() for pp in ppb.build_peptides(chain)),Seq("",IUPAC.protein))
            seq_spec = None #TODO: implement
        elif method == "order":
            seq = []
            seq_spec = []
            for res in chain.get_residues():
                seq.append(three_to_one(res.get_resname()))
                ## from Bio docs, res.get_full_id() returns: ("1abc", 0, "A", (" ", 10, "A"))
                fid = res.get_full_id()
                seq_spec.append(pdb_seq_spec(
                    chain = fid[-2].strip(),
                    resn = res.get_resname(),
                    resi = fid[-1][-2],
                    ins = fid[-1][-1].strip()
                    ))


            seq = Seq("".join(seq),IUPAC.protein)
        else:
            raise ValueError("Unknown method: {}".format(method))

        seq_chains.append(dict(
            id_chain=id_chain, 
            seq_rec=SeqRecord(seq,id="{}_{}".format(id,id_chain),description=""),
            seq_spec=seq_spec
            ))
        chains_map = dict(((x["id_chain"],x) for x in seq_chains))
    return pdb_seqs(id=id,chains=seq_chains,chains_map=chains_map)


def pdb_sequence_save(pdb_file,id=None,fasta_file=None,strip_sfx="."):

    if id is None:
        id = util.make_id_from_file_name(pdb_file)
    if fasta_file is None:
        fasta_file = id+".fasta"
    print(fasta_file)
    with open(fasta_file,"w") as out_fa:
        seq_res = pdb_sequence(pdb_file=pdb_file,id=id)
        for chain in seq_res["chains"]:
            SeqIO.write(chain["seq_rec"],out_fa,"fasta")


def fasta_to_df(seq_file,seq_format_out="str",alphabet=None):
    from Bio import SeqIO
    assert seq_format_out in ("seq","str"), "Unknown seq_format_out: {}".format(seq_format_out)
    with util.open_gzip_text23(seq_file,"rt") as inp:
        recs = [dict(id=rec.id,description=rec.description,seq=str(rec.seq) \
            if seq_format_out == "str" else rec.seq) \
            for rec in SeqIO.parse(inp, "fasta",alphabet=alphabet)]
        return pd.DataFrame(recs) 


def _main():
    import argh
    argh.dispatch_commands([pdb_sequence_save])

