from __future__ import unicode_literals
from __future__ import print_function
from __future__ import division
from __future__ import absolute_import
from future import standard_library
standard_library.install_aliases()
from builtins import *

import helpers

import pytest
from subprocess import check_call, check_output
import os
import shlex
import shutil

from MGT.Logging import *

pytestmark = pytest.mark.usefixtures("goto_cleandir_test","get_test_data_dir")

test_data = "test_data"
test_data_cwl_toil = os.path.join(test_data,"cwl/toil")


@helpers.skip_no_conda
def test_init_dep_resolver():
    cmd = "python -m MICGENT.cwl_runner add-dep-resolver --pkg-name ngs-mstb --conda-env ariba --mode w my_dep_resolver.yaml"
    check_call(shlex.split(cmd))
    assert os.path.isfile("my_dep_resolver.yaml")
    assert os.path.isfile("deps/ngs-mstb/1.0/env.sh")
    assert helpers.grep(r"ariba","deps/ngs-mstb/1.0/env.sh")

@helpers.skip_no_conda
def test_append_dep_resolver():
    cmd = "python -m MICGENT.cwl_runner add-dep-resolver --pkg-name ngs-mstb --conda-env ariba --mode w my_dep_resolver.yaml"
    check_call(shlex.split(cmd))
    cmd = "python -m MICGENT.cwl_runner add-dep-resolver --pkg-name ngs-mstb-py2 --conda-env ariba my_dep_resolver.yaml"
    check_call(shlex.split(cmd))
    assert os.path.isfile("my_dep_resolver.yaml")
    assert os.path.isfile("deps/ngs-mstb-py2/1.0/env.sh")
    assert helpers.grep(r"ariba","deps/ngs-mstb-py2/1.0/env.sh")
    ## previous env is still there?
    assert os.path.isfile("deps/ngs-mstb/1.0/env.sh")
    assert helpers.grep(r"ariba","deps/ngs-mstb/1.0/env.sh")


@helpers.skip_no_conda_env_toil
def test_run_toil():
    check_call("cp {}/* ./".format(test_data_cwl_toil),shell=True)
    cmd = "python -m MICGENT.cwl_runner run-toil --logLevel DEBUG --runner-use-conda --runner-conda-env toil sorttool.cwl revsort-job.json"
    check_call(shlex.split(cmd))
    assert os.path.isfile("out/output.txt")
    assert helpers.grep(r"whenever it is a damp","out/output.txt")
