from __future__ import unicode_literals
from __future__ import print_function
from __future__ import division
from __future__ import absolute_import
from future import standard_library
standard_library.install_aliases()
from builtins import *

from MICGENT import arg_parsing
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
test_data_util = os.path.join(test_data,"util")

pjoin = os.path.join

@argh.arg("--config",type=yaml_util.load_yaml)
def my_g(x,y="a",config={},out="my_g_out.yaml"):
    ret = dict(x=x,y=y,config=config,out=out)
    return ret


@argh.arg("--config",type=yaml_util.load_yaml)
def my_g_exp_1(x=15,y="b",config={},out="my_g_out.yaml"):
    ret = dict(x=x,y=y,config=config,out=out)
    return ret

@argh.arg("--config",type=yaml_util.load_yaml)
def my_g_exp_2(x=17,y="P",config={},out="my_g_out.yaml"):
    ret = dict(x=x,y=y,config=config,out=out)
    return ret

def test_arg_parsing_cli():
    cli_script = pjoin(test_data_util,"arg_parsing_test_cli.py")
    conf = pjoin(test_data_util,"arg_parsing_test_cli.yaml")
    conf_out_exp = pjoin(test_data_util,"arg_parsing_test_cli_out_exp.yaml")
    conf_out = "arg_parsing_test_cli_out.yaml"
    cmd = "python {cli_script} --config {conf} my-g --config {conf} 17".format(**locals())
    check_call(shlex.split(cmd))
    assert os.path.isfile(conf_out)
    conf_out_o = yaml_util.load_yaml(conf_out)
    conf_out_exp_o = yaml_util.load_yaml(conf_out_exp)
    assert conf_out_o == conf_out_exp_o

def test_arg_parsing_api():
    my_g_par = arg_parsing.update_sig(my_g, x=15, y="b")

    log.info("Updated signature: {}".format(argh.utils.get_arg_spec(my_g_par)))

    assert my_g_par(x=15, y="b") == my_g_exp_1()

    conf = dict(x=17, y="P")

    my_g_par = arg_parsing.update_sig(my_g, **conf)

    log.info("Updated signature: {}".format(argh.utils.get_arg_spec(my_g_par)))

    assert my_g_par(x=17, y="P") == my_g_exp_2()
