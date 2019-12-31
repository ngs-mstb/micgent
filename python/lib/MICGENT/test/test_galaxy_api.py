from __future__ import unicode_literals
from __future__ import print_function
from __future__ import division
from __future__ import absolute_import
from future import standard_library
standard_library.install_aliases()
from builtins import *

from MICGENT.galaxy import api as gal_api
from MICGENT import yaml_util
import argh

import helpers

import pytest
from subprocess import check_call, check_output
import os
import shlex
import shutil

from MGT.Logging import *

pytestmark = pytest.mark.usefixtures("goto_cleandir_test","get_test_data_dir")

test_data = "test_data"
test_data_galaxy = os.path.join(test_data,"galaxy")

pjoin = os.path.join

@pytest.mark.skip(reason="This test needs live Galaxy instance and respective arguments. Not implemented yet.")
def test_simplify_rerun_json():
    test_data_galaxy = os.path.abspath(globals()["test_data_galaxy"])
    job_yaml_exp = pjoin(test_data_galaxy,"job.yaml")
    job_user_yaml_exp = pjoin(test_data_galaxy,"job_user.yaml")
    job_exp = yaml_util.load_yaml(job_yaml_exp)
    with helpers.mkchdir("simplify_rerun_json"):
        ## job_res will be a totaly independent identical copy of job_exp;
        ## job_res will be passed as input argument to simplify_rerun_json,
        ## and we expect it will not be modified as a side-effect
        job_res = yaml_util.load_yaml(job_yaml_exp)
        job_user_exp = yaml_util.load_yaml(job_user_yaml_exp)
        job_user_res = gal_api.simplify_rerun_json(job_res)
        yaml_util.dump_yaml(job_user_res,"job_user.yaml")
        if job_exp != job_res:
            print(yaml_util.dumps_yaml(yaml_util.diff_yaml(job_exp,job_res)))
            raise AssertionError("Job YAML object does not match expected value")
        if job_user_exp != job_user_res:
            print(yaml_util.dumps_yaml(yaml_util.diff_yaml(job_user_exp,job_user_res)))
            raise AssertionError("Job YAML object does not match expected value")
