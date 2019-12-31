"""Methods for processing annotations"""
from __future__ import unicode_literals
from __future__ import print_function
from __future__ import division
from __future__ import absolute_import
from future import standard_library
standard_library.install_aliases()
from builtins import *
from MGT.Logging import *

def select_top_overlapping_intervals(begin,end,closed_end=False):
    """Given a collection of overlapping intervals, select one "best" for each overlapping subset.
    For example, intervals can be BLASTN hits from reference DB to the contig.

    Intervals should be provided as begin and end iterables, **sorted in decreasing** order
    by "interval strength" that will define the preference when selecting
    among overlapping intervals.

    result - this is a generator that yields indices of the selected records in the input.

    How selection is currently done: if A,B and C are intervals,
    shown here is decreasing order of preference, and B overlaps both
    A and C, but A does not overlap C, then both A and C will be returned
    as selections. Moving in the initial order of intervals, the algorithm
    discards from future consideration all intervals that overlap with the
    currently considered intervals, stores the index of the current interval,
    and moves to the next remaining interval.
    """
    ##If other variations of pruning the overlapping intervals are needed,
    ##we can use any of the Python graph libraries such as APGL, build a
    ##graph from overlaps and use methods like findConnectedComponents() etc.
    from intervaltree import Interval, IntervalTree
    intervals = [Interval(*iv) for iv in zip(begin,
                                                 end if not closed_end else ((_+1) for _ in end))]
    tree = IntervalTree(intervals)
    ##TODO: what is the cost of removing tree nodes? Maybe it is cheaper
    ##to tag nodes overlapping the current one is a separate array.
    for iiv,iv in enumerate(intervals):
        ## docs says tree membership check for Interval object is O(1)
        if iv in tree:
            yield iiv
            try:
                tree.remove_overlap(iv.begin,iv.end)
            except KeyError as msg:
                log.warning("KeyError when removing existing node from IntervalTree: {}. \
                        This must be a bug to fix in IntervalTree code. Ignoring this error for now.".\
                        format(msg))

