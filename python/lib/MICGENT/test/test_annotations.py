from __future__ import unicode_literals
from __future__ import print_function
from __future__ import division
from __future__ import absolute_import
from future import standard_library
standard_library.install_aliases()
from builtins import *
from MGT.Logging import *

def test_select_top_overlapping_intervals():
    import MICGENT.annotations
    import pandas as pd
    ivs = pd.DataFrame.from_records([(1, 10),
                                     (10, 20),
                                     (11, 40),
                                     (100,110)],
                              columns=["begin", "end"])

    sel_ind = list(MICGENT.annotations.select_top_overlapping_intervals(ivs.begin,ivs.end,closed_end=True))
    log.info(sel_ind)
    assert sel_ind == [0, 2, 3]

    sel_ind = list(MICGENT.annotations.select_top_overlapping_intervals(ivs.begin, ivs.end, closed_end=False))
    log.info(sel_ind)
    assert sel_ind == [0, 1, 3]
