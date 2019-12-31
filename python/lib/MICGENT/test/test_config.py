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

jsonnet = pytest.importorskip("jsonnet")

pytestmark = pytest.mark.usefixtures("goto_cleandir_test","get_test_data_dir")

test_data = "test_data"

from MGT.Logging import *

def test_load_config_from_package_vars_default():
    import MICGENT.config
    conf = MICGENT.config.load_config(pkg="MICGENT",
                                        vars_default=dict(data="some_data"))
    log.debug(conf)

def test_load_config_from_package_environ_vars():
    import MICGENT.config
    os.environ["MICGENT_DATA"] = "some_data"
    conf = MICGENT.config.load_config(pkg="MICGENT")
    log.debug(conf)


def test_load_config_from_environ():
    import MICGENT.config
    os.environ["MICGENT_CONFIG"] = pjoin(test_data,"sample_config.jsonnet")
    os.environ["MICGENT_DATA"] = "some_data"

    seq_db = MICGENT.config.load_config(prefix_env="MICGENT")["data"]["cge"]["plasmid_db"]["seq"]
    log.debug(seq_db)
    assert seq_db == "some_data/cge/plasmid_finder/2016-01-04/plasmidfinder/plasmid_database.fsa"
