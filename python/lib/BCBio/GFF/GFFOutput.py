"""Output Biopython SeqRecords and SeqFeatures to GFF3 format.

The target format is GFF3, the current GFF standard:
    http://www.sequenceontology.org/gff3.shtml
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
from copy import copy
import textwrap
import re

from Bio import SeqIO

class _IdHandler(object):
    """Generate IDs for GFF3 Parent/Child relationships where they don't exist.
    """
    def __init__(self):
        self._prefix = "biopygen"
        self._counter = 1
        self._seen_ids = set()

    def _generate_id(self, quals):
        """Generate a unique ID not present in our existing IDs.
        """
        gen_id = self._get_standard_id(quals)
        if gen_id and gen_id not in self._seen_ids:
            return gen_id
        if gen_id is None:
            prefix = self._prefix
        else:
            prefix = gen_id
        while 1:
            gen_id = "%s.%s" % (prefix, self._counter)
            if gen_id not in self._seen_ids:
                break
            self._counter += 1
        return gen_id

    def _get_standard_id(self, quals):
        """Retrieve standardized IDs from other sources like NCBI GenBank.

        This tries to find IDs from known key/values when stored differently
        than GFF3 specifications.
        """
        cur_id = None
        possible_keys = ["transcript_id", "protein_id", "Name"]
        for test_key in possible_keys:
            if test_key in quals:
                cur_id = quals[test_key]
                if isinstance(cur_id, tuple) or isinstance(cur_id, list):
                    if cur_id:
                        cur_id = cur_id[0]
                    else:
                        cur_id = None
        return cur_id

    def update_quals(self, quals, has_children):
        """Update a set of qualifiers, adding an ID if necessary.
        """
        cur_id = quals.get("ID", None)
        # if we have an ID, record it
        if cur_id:
            if not isinstance(cur_id, list) and not isinstance(cur_id, tuple):
                cur_id = [cur_id]
            for add_id in cur_id:
                self._seen_ids.add(add_id)
        # if we need one and don't have it, create a new one
        elif has_children:
            new_id = self._generate_id(quals)
            self._seen_ids.add(new_id)
            quals["ID"] = [new_id]
        return quals


_non_quote = " a-zA-Z0-9.:\^*$@!+_?-|"
# official (from gmod), but space is only not allowed in ID: _non_quote="[a-zA-Z0-9.:^*$@!+_?-|]
def quote_gff_val(x):
    import six
    return six.moves.urllib.parse.quote(str(x), safe=_non_quote)

class GFF3Writer(object):
    """Write GFF3 files starting with standard Biopython objects.
    """
    def __init__(self):
        pass

    def write(self, recs, out_handle, include_fasta=False, out_handle_fasta=None, rec_id_field="id"):
        """Write the provided records to the given handle in GFF3 format.
        """
        id_handler = _IdHandler()
        self._write_header(out_handle)
        fasta_recs = []
        try:
            recs = iter(recs)
        except TypeError:
            recs = [recs]
        for rec in recs:
            self.last_by_type = {}
            if not rec.id:
                rec.id = rec.name
            rec.id = getattr(rec,rec_id_field)
            self._write_rec(rec, out_handle)
            self._write_annotations(rec.annotations, rec.id, len(rec.seq), out_handle)
            for sf in rec.features:
                sf = self._clean_feature(sf)
                self.last_by_type[sf.type] = sf
                id_handler = self._write_feature(sf, rec.id, out_handle,
                        id_handler)
            if include_fasta and len(rec.seq) > 0:
                fasta_recs.append(rec)
        if len(fasta_recs) > 0:
            fa_handle = out_handle
            fasta_destination = "genbank"
            if out_handle_fasta:
                fa_handle = out_handle_fasta
                fasta_destination = "fasta"
            self._write_fasta(fasta_recs, fa_handle, destination=fasta_destination)

    def _clean_feature(self, feature):
        quals = {}
        for key, val in list(feature.qualifiers.items()):
            if not isinstance(val, (list, tuple)):
                val = [val]
            val = [str(x) for x in val]
            quals[key] = val
        feature.qualifiers = quals
        return feature

    def _write_rec(self, rec, out_handle):
        # if we have a SeqRecord, write out optional directive
        if len(rec.seq) > 0:
            out_handle.write("##sequence-region %s 1 %s\n" % (rec.id, len(rec.seq)))

    def _get_phase(self, feature):
        if "phase" in feature.qualifiers:
            phase = feature.qualifiers["phase"][0]
        elif feature.type == "CDS":
            phase = int(feature.qualifiers.get("codon_start", [1])[0]) - 1
        else:
            phase = "."
        return str(phase)

    def _write_feature(self, feature, rec_id, out_handle, id_handler,
            parent_id=None):
        """Write a feature with location information.
        """
        def _select_first_qual_non_empty(quals,keys):
            for key in keys:
                if quals.get(key,None):
                    return key,quals[key]
            return None,list()
        if feature.strand == 1:
            strand = '+'
        elif feature.strand == -1:
            strand = '-'
        else:
            strand = '.'
        # remove any standard features from the qualifiers
        quals = feature.qualifiers.copy()
        if 'translation' in quals:
            quals['translation'] = [ re.sub(r'\s','',x) for x in quals['translation'] ]
        for std_qual in ["source", "score", "phase"]:
            if std_qual in quals and len(quals[std_qual]) == 1:
                del quals[std_qual]
        # add a link to a parent identifier if it exists
        if parent_id:
            if "Parent" not in quals:
                quals["Parent"] = []
            quals["Parent"].append(parent_id)
        if not quals.get("Name",[]):
            quals["Name"] = _select_first_qual_non_empty(quals,("product","gene","organism"))[1]
        quals = id_handler.update_quals(quals, True)
        feature.qualifiers["ID"] = quals["ID"]

        if feature.type:
            ftype = feature.type
        else:
            ftype = "sequence_feature"
        if ftype == "mat_peptide":
            last_CDS = self.last_by_type.get("CDS",None)
            if last_CDS:
                gene_CDS = last_CDS.qualifiers.get("gene",None)
                gene_this = quals.get("gene",None)
                if gene_this is not None and gene_this == gene_CDS:
                    #we cannot use the correct Derives_from qual because
                    #JBrowse will not show mat_peptides by instead create
                    #ugly text box attribute Derived Features for the parent CDS
                    quals["Parent_CDS"] = last_CDS.qualifiers["ID"]
        parts_loc_ind = (3,4)
        parts = [str(rec_id),
                 feature.qualifiers.get("source", ["feature"])[0],
                 ftype,
                 str(feature.location.nofuzzy_start + 1), # 1-based indexing
                 str(feature.location.nofuzzy_end),
                 feature.qualifiers.get("score", ["."])[0],
                 strand,
                 self._get_phase(feature),
                 self._format_keyvals(quals)]
        out_handle.write("\t".join(parts) + "\n")

        if len(feature.location.parts) >= 2 or feature.type == "CDS":
            for sub_location in feature.location.parts:
                sub_feature = copy(feature)
                sub_par_id = quals["ID"][0]
                sub_quals = copy(quals)
                sub_feature.qualifiers = sub_quals
                del sub_quals["ID"]
                for key in list(sub_feature.qualifiers.keys()):
                    if key not in ("Name","ID"):
                        del sub_feature.qualifiers[key]
                sub_feature.location = sub_location
                if sub_feature.type == "CDS":
                    sub_feature.type = "exon"
                    if "translation" in sub_feature.qualifiers:
                        del sub_feature.qualifiers["translation"]

                id_handler = self._write_feature(sub_feature, rec_id, out_handle,
                        id_handler, sub_par_id)
        return id_handler

    def _format_keyvals(self, keyvals):
        format_kvs = []
        for key, values in list(keyvals.items()):
            key = key.strip()
            format_vals = []
            if not isinstance(values, list) or isinstance(values, tuple):
                values = [values]
            for val in values:
                s_val = str(val).strip()
                #s_val = textwrap.fill(s_val,60)
                val = quote_gff_val(s_val)
                if ((key and val) and val not in format_vals):
                    format_vals.append(val)
            #JBrowse does not show attributes w/o values, so we add a space
            if len(format_vals) == 0:
                format_vals.append(quote_gff_val(" "))
            format_kvs.append("%s=%s" % (key, ",".join(format_vals)))
        return ";".join(format_kvs)

    def _write_annotations(self, anns, rec_id, size, out_handle):
        """Add annotations which refer to an entire sequence.
        """
        format_anns = self._format_keyvals(anns)
        if format_anns:
            parts = [rec_id, "annotation", "remark", "1", str(size if size > 1 else 1),
                     ".", ".", ".", format_anns]
            out_handle.write("\t".join(parts) + "\n")

    def _write_header(self, out_handle):
        """Write out standard header directives.
        """
        out_handle.write("##gff-version 3\n")

    def _write_fasta(self, recs, out_handle, destination="genbank"):
        """Write sequence records using the ##FASTA directive.
        """
        if destination == "genbank":
            out_handle.write("##FASTA\n")
        SeqIO.write(recs, out_handle, "fasta")

def write(recs, out_handle, include_fasta=False, out_handle_fasta=None, rec_id_field="id"):
    """High level interface to write GFF3 files from SeqRecords and SeqFeatures.

    If include_fasta is True, the GFF3 file will include sequence information
    using the ##FASTA directive.
    """
    writer = GFF3Writer()
    return writer.write(recs, out_handle, include_fasta=include_fasta,
                        out_handle_fasta=out_handle_fasta, rec_id_field=rec_id_field)
