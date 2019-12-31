from __future__ import unicode_literals
from __future__ import print_function
from __future__ import division
from __future__ import absolute_import
from future import standard_library
standard_library.install_aliases()
from builtins import *
import pytest
from subprocess import check_call
import os, glob
from os.path import join as pjoin

pytestmark = pytest.mark.usefixtures("goto_cleandir_test","get_test_data_dir")


test_data = "test_data"
test_data_plasmid = pjoin(test_data,"plasmid_finder")
test_data_res = pjoin(test_data,"res_finder")

from MGT.Logging import *

@pytest.mark.xfail(reason="Need to update the expected output")
def test_annotate_plasmids():
    import MICGENT.finders_cge
    import pandas as pd
    bl_df = MICGENT.finders_cge.annotate_plasmids(seq=pjoin(test_data_plasmid,"test.query.fsa"),
                                                     seq_db=pjoin(test_data_plasmid,"test.db.fsa"),
                                                     db_cov_min = 0.0)
    log.info(bl_df)
    assert bl_df[["query_id","sbjct_id"]].equals(
        pd.DataFrame.from_records([("Q1","D2"),
                                   ("Q1","D1"),
                                   ("Q1","D1"),
                                   ("Q2","D2"),
                                   ("Q2","D2"),
                                   ("Q3","D3"),
                                   ("Q3","D3")],
                                  columns=["query_id","sbjct_id"])
    )

    bl_df = MICGENT.finders_cge.annotate_plasmids(seq=pjoin(test_data_plasmid,"test.query.fsa"),
                                                     seq_db=pjoin(test_data_plasmid,"test.db.fsa"),
                                                     db_cov_min = 0.6)
    log.info(bl_df)
    assert bl_df[["query_id","sbjct_id"]].equals(
        pd.DataFrame.from_records([("Q2","D2")],
                                  columns=["query_id","sbjct_id"])
    )

def test_annotate_resistance_genes():
    import MICGENT.finders_cge
    import pandas as pd
    bl_df = MICGENT.finders_cge.annotate_resistance_genes(seq=pjoin(test_data_res,"test.query.fsa"),
                                                     seq_db=pjoin(test_data_res,"test.db.?.fsa"))
    log.info(bl_df)
    assert bl_df[["query_id","sbjct_id","res_class"]].equals(
        pd.DataFrame.from_records([("blaTMB-1_1_FR771847","blaTMB-1_1_FR771847","blaTMB-1_1"),
                                   ("blaTMB-2_1_AB758277","blaTMB-2_1_AB758277","blaTMB-2_1")],
                                  columns=["query_id","sbjct_id","res_class"])
    )

@pytest.mark.micgent_data
def test_annotate_resistance_genes_full_db():
    _jsonnet = pytest.importorskip("_jsonnet")
    import MICGENT.finders_cge
    bl_df = MICGENT.finders_cge.annotate_resistance_genes(seq=pjoin(test_data_res,"test.query.fsa"))
    log.info(bl_df)

def test_select_top_overlapping_hits_cge():
    import MICGENT.finders_cge
    import pandas as pd
    hits = [
        ["ind","query_id", "length_score", "pct_ident", "sbjct_id", "query_start", "query_end","locus"],
        [0, "Q1", 1, 90,  "D2", 1,  10, "Locus1"],
        [1, "Q1", 2, 100, "D3", 10, 20, "Locus1"],
        [2, "Q1", 2, 95,  "D4", 11, 40, "Locus1"],
        [3, "Q1", 3, 80,  "D5", 100,110, "Locus1"],
        [4, "Q2", 1, 90,  "D2", 1,  10, "Locus1"],
        [5, "Q2", 2, 100, "D3", 10, 20, "Locus1"],
        [6, "Q2", 2, 95,  "D4", 11, 40, "Locus1"],
        [7, "Q3", 3, 80,  "D5", 100,110, "Locus1"],
    ]
    hits = pd.DataFrame.from_records(hits[1:],
                              columns=hits[0])
    hits_sel = MICGENT.finders_cge.select_top_overlapping_hits_cge(hits.sample(frac=1).reset_index(drop=True),
                                                                   closed_end=True)
    log.info(hits_sel)
    hits_expect = hits.iloc[[0, 2, 3, 4, 6, 7],:]
    hits_expect.reset_index(drop=True,inplace=True)
    assert hits_sel.equals(hits_expect)
