"""Top level of GFF parsing providing shortcuts for useful classes.
"""
from __future__ import absolute_import
from __future__ import unicode_literals
from __future__ import print_function
from __future__ import division
from future import standard_library
standard_library.install_aliases()
from builtins import *
from .GFFParser import GFFParser, DiscoGFFParser, GFFExaminer, parse, parse_simple
from .GFFOutput import GFF3Writer, write
