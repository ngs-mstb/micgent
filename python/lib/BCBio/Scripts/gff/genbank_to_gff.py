#!/usr/bin/env python
"""Convert a GenBank file into GFF format.

Usage:
    genbank_to_gff.py <genbank_file>
"""
from __future__ import unicode_literals
from __future__ import print_function
from __future__ import division
from __future__ import absolute_import
from future import standard_library
standard_library.install_aliases()
from builtins import *
import sys
import os

from Bio import SeqIO
from Bio import Seq

from BCBio import GFF

def main(gb_file,include_fasta=None):
    out_file = "%s.gff" % os.path.splitext(gb_file)[0]
    inc_fasta = False
    if include_fasta is not None:
        if include_fasta.lower() in ("true","yes","1"):
            inc_fasta = True
        
    with open(out_file, "w") as out_handle:
        GFF.write(SeqIO.parse(gb_file, "genbank"), out_handle, inc_fasta)

if __name__ == "__main__":
    main(*sys.argv[1:])
