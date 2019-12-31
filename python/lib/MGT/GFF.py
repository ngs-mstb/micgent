"""Support for writing GFF3 files.
The format has a short description here:
http://gmod.org/wiki/GFF3
A formal definition:
http://www.sequenceontology.org/gff3.shtml
And the formal definition of attributes("ontology") in SOFA reference on this page:
http://www.sequenceontology.org/resources/intro.html
"""
from __future__ import unicode_literals
from __future__ import print_function
from __future__ import division
from __future__ import absolute_import
from future import standard_library
standard_library.install_aliases()
from builtins import str
from builtins import *
from builtins import object
import re
import six

class GFF3Line(object):
    """Abstract class for GFF3 line of output"""
    def write(self,out):
        out.write(str(self))

class GFF3Pragma(GFF3Line):
    """GFF3 pragma abstract class"""

    def __str__(self):
        v_s = self.get_value_str()
        n_s = "##{}".format(self.name)
        if v_s is not None:
            return "{}\t{}\n".format(n_s,v_s)
        else:
            return "{}\n".format(n_s)

    def get_value_str(self):
        return ""


class GFF3Header(GFF3Pragma):
    """GFF3 file header pragma"""

    name = "gff-version"

    def get_value_str(self):
        return "3.2.1"

class GFF3Fasta(GFF3Pragma):
    """GFF3 file FASTA pragma"""

    name = "FASTA"

    def get_value_str(self):
        return None


class GFF3SeqRegion(GFF3Pragma):
    """GFF3 file sequence region pragma"""

    name = "sequence-region"

    def __init__(self,seqid,start,end):

        self.seqid = seqid
        self.start = start
        self.end = end

    def get_value_str(self):
        return "{}\t{}\t{}".format(self.seqid,self.start,self.end)



class GFF3Attributes(dict):
    """Represents GFF3 attributes."""

    _non_quote=" a-zA-Z0-9.:\^*$@!+_?-|"
    #official (from gmod), but space is only not allowed in ID: _non_quote="[a-zA-Z0-9.:^*$@!+_?-|]

    _re_quote=re.compile("[^%s]" % _non_quote)

    @classmethod
    def quote_val_url(klass,val):
        return six.moves.urllib.parse.quote(str(val),safe=klass._non_quote)

    @classmethod
    def quote_val_re(klass,val):
        return re.sub(klass._re_quote," ",str(val))

    def __init__(self,_quote_method="url",**kw):
        dict.__init__(self,**kw)
        ## set quote_method to select how characters will be quoted
        ## when writing GFF3 file. Although the standard tells to use URL quoting,
        ## some tools (e.g. genometools) do not unquote them when creating diagrams,
        ## which lead to ugly looking text.
        ## Posiible values are:
        ## re  - replace every non-allowed symbol with space
        ## url [default] - use url quoting
        self._quote_method = _quote_method

    def __str__(self):
        if self._quote_method == "url":
            quote_val = self.quote_val_url
        else:
            quote_val = self.quote_val_re
        return ';'.join(( "%s=%s" % (tag, quote_val(value) if isinstance(value,six.string_types) or \
                not hasattr(value,"__len__") else ','.join( (quote_val(x) for x in value) )) for \
                (tag,value) in sorted(self.items())))

GFF3At = GFF3Attributes

class GFF3Record(GFF3Line):
    """One feature line in GFF3 file."""


    def __init__(self,seqid='.',source='.',type='.',start=0,end=1,
            score='.',strand='.',phase='.',attribs=None,**kw):
        """Constructor.
        Attributes can be passed as GFF3Attributes or dict object through attribs argument, or
        as additional keyword arguments, or both.
        @param start Start of feature (zero-based, will be converted to GFF unit-based on output)
        @end end of feature (zero-based, will be converted to GFF unit-based on output)"""
        self.seqid = seqid
        self.source = source
        self.type = type
        self.start = start
        self.end = end
        self.score = score
        if isinstance(strand,(int,float)):
            if strand == 0:
                strand = "."
            elif strand < 0:
                strand = "-"
            else:
                strand = "+"
        self.strand = strand
        if type == "CDS" and phase == ".":
            ##for CDS, phase should always be defined
            phase = 0
        self.phase = phase

        if attribs is None:
            self.attribs = GFF3Attributes(**kw)
        else:
            self.attribs = GFF3Attributes(**attribs)
            self.attribs.update(kw)

    def __str__(self):
        return "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (self.seqid,self.source,self.type,self.start+1,self.end,
                self.score,self.strand,self.phase,self.attribs)

    def copy(self):
        from copy import copy
        c = copy(self)
        c.attribs = copy(self.attribs)
        return c

    def clone(self):
        """Creates a copy of the current object with empty attributes list"""
        c = self.copy()
        c.attribs = GFF3Attributes()

    def fromSeqFeat(self,feat):
        """Pull data from Bio.SeqFeature instance.
        This mutates the current object"""
        self.type = feat.type
        self.start = feat.location.nofuzzy_start
        self.end = feat.location.nofuzzy_end
        self.strand = '+' if feat.strand > 0 else '-' if feat.strand < 0 else '.'
        # no attributes should be automatically carried forward from the previous values:
        self.attribs = GFF3Attributes()
        ats = self.attribs
        quals = feat.qualifiers
        ret = [ self ]
        if "id" in quals:
            ats["ID"] = quals["id"][0]
        if self.type == "CDS":
            self.phase = 0
            ats["Name"] = quals["product"][0]
            ats["ID"] = quals["protein_id"]
        elif self.type == "gene":
            ats["Name"] = quals["locus_tag"][0]
            ats["ID"] = ats["Name"]
        else:
            ats["Name"] = quals.get("note",quals.get("id",(self.type,)))[0]
        if quals.get("rpt_type",("",))[0] == "CRISPR":
            self.type = "CRISPR"
            self.strand = '.'
            ats["ID"] = quals["note"]
            if "rpt_unit_range" in quals:
                for srange in quals["rpt_unit_range"]:
                    start,end = srange.split("..")
                    start = int(start)-1
                    end = int(end)
                    rec = self.copy()
                    rec.type = "repeat_unit"
                    rec.start = start
                    rec.end = end
                    rec_ats = rec.attribs
                    del rec_ats["Name"], rec_ats["ID"]
                    rec_ats["Parent"] = ats["ID"]
                    ret.append(rec)
        #if ats["ID"] == "ACJ74823.1":
        #    pdb.set_trace()
        return ret


