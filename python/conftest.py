from __future__ import unicode_literals
from __future__ import print_function
from __future__ import division
from __future__ import absolute_import

from future import standard_library
standard_library.install_aliases()
from builtins import *
import pytest
import tempfile, os, shutil, importlib

from MICGENT import resources

_pkg_name = "MICGENT"

@pytest.fixture(scope="session")
def get_test_data_dir():
    importlib.import_module(".resources",_pkg_name)
    inp_test_data = resources.get_pkg_test_data_dir(_pkg_name)
    shutil.copytree(inp_test_data,os.path.join(os.getcwd(),"test_data"))

@pytest.fixture(scope="session")
def goto_cleandir_test():
    root_test_dir = "test_run"
    if not os.path.exists(root_test_dir):
        os.makedirs(root_test_dir)
    newpath = tempfile.mkdtemp(prefix="run.test.",suffix=".tmp",dir=root_test_dir)
    os.chdir(newpath)

## Note that you need to use = between the long option name and value, otherwise you
## will get unrecongnized argument errors when using both options (it interprets
## space-separated values as extra options)
## `pytest --micgent-data=~/work/micgent_db --run-slow`
## And `conftest.py` has to be located in the directory where you are running `pytest`
## command - apparently pytest does not discover deeply placed `conftest.py` files
## before trying to interpret the command arguments.
def pytest_addoption(parser):
    parser.addoption('--micgent-data', action='store', default=None, help="Path to MICGENT DB")
    parser.addoption('--large-test-data', action='store', default=None, help="Path to MICGENT large test data",
                     type=os.path.expanduser)
    parser.addoption('--huge-test-data', action='store', default=None, help="Path to MICGENT huge test data",
                     type=os.path.expanduser)
    parser.addoption('--conda-env-ngs-mstb', action='store', default=None, help="Name of existing Conda environment "
                                                                                "for ngs-mstb Gene Extractor CWL package dependency")
    parser.addoption("--run-slow", action="store_true",
                     default=False, help="run slow tests")
    parser.addoption('--extra-config', action='store', default=None, help="Path to the extra config YAML in the format of MICGENT.arg_parsing module",
                     type=lambda x: os.path.abspath(os.path.expanduser(x)))

def pytest_collection_modifyitems(config, items):
    micgent_data = os.path.expanduser(config.getoption("--micgent-data"))
    if micgent_data:
        # --micgent-data given in cli: do not skip tests that need it
        # set env var for finders_cge to find the DB
        os.environ["MICGENT_DATA"] = micgent_data
    else:
        skip_no_micgent_data = pytest.mark.skip(reason="need --micgent-data option to run")
        for item in items:
            if "micgent_data" in item.keywords:
                item.add_marker(skip_no_micgent_data)
    if not config.getoption("--large-test-data"):
        skip_large_test_data = pytest.mark.skip(reason="need --large-test-data option to run")
        for item in items:
            if "large_test_data" in item.keywords:
                item.add_marker(skip_large_test_data)
    if not config.getoption("--huge-test-data"):
        skip_huge_test_data = pytest.mark.skip(reason="need --huge-test-data option to run")
        for item in items:
            if "huge_test_data" in item.keywords:
                item.add_marker(skip_huge_test_data)
    if not config.getoption("--run-slow"):
        skip_slow = pytest.mark.skip(reason="need --run-slow option to run")
        for item in items:
            if "slow" in item.keywords:
                item.add_marker(skip_slow)
    if not config.getoption("--conda-env-ngs-mstb"):
        skip_conda_env_ngs_mstb = pytest.mark.skip(reason="need --conda-env-ngs-mstb option to run")
        for item in items:
            if "conda_env_ngs_mstb" in item.keywords:
                item.add_marker(skip_conda_env_ngs_mstb)
