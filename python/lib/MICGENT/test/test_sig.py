from __future__ import unicode_literals
from __future__ import print_function
from __future__ import division
from __future__ import absolute_import
from future import standard_library
standard_library.install_aliases()
from builtins import *

from MICGENT import sig
from MICGENT import yaml_util

import helpers

import pytest
import subprocess
import os
import shlex
import shutil

from MGT.Logging import *

pytestmark = pytest.mark.usefixtures("goto_cleandir_test","get_test_data_dir")

test_data = "test_data"
test_data_util = os.path.join(test_data,"util")

pjoin = os.path.join

def test_file_sig_cmp_msg():
    """Test config file interface and error message generation in signature comparisons"""
    test_fn = "test.txt"
    sig_conf_fn = "sig.yaml"
    test_sig_fn = "test.sig"
    key_other = "65757ffsswgw8g"
    with open(test_fn,"w") as out:
        out.write("FFUGJGJHGJGYNVNVvhgjgjyg76868uyu\n")
    yaml_util.dump_yaml(dict(file_sig_cmp_msg=dict(key=key_other)),sig_conf_fn)
    test_sig = sig.file_sig(test_fn,out_file=test_sig_fn)
    sig.file_sig_cmp_msg(sig_file=test_sig_fn,inp_file=test_fn)
    with pytest.raises(SystemExit):
        sig.file_sig_cmp_msg(sig_file=test_sig_fn,inp_file=test_fn,key=key_other)
    with subprocess.Popen(["python","-m","MICGENT.sig","--config",sig_conf_fn,"file-sig-cmp-msg",
        "--msg","Sig mismatch",test_sig_fn,test_fn],stderr=subprocess.PIPE,universal_newlines=True) as p:
        p.wait() # works around bug in Python 3.6 context manager - can get p.returncode None
        assert p.returncode != 0, "Expected a return code indicating caught error"
        assert p.stderr.read().strip() == "Sig mismatch", "Error message is not what was expected"
    sig.file_sig(test_fn,key=key_other,out_file=test_sig_fn)
    with subprocess.Popen(["python","-m","MICGENT.sig","--config",sig_conf_fn,"file-sig-cmp-msg",
        "--msg","Sig mismatch",test_sig_fn,test_fn],stderr=subprocess.PIPE,universal_newlines=True) as p:
        p.wait()  # works around bug in Python 3.6 context manager - can get p.returncode None
        assert p.returncode == 0, "Expected a return code indicating caught error"        
        assert p.stderr.read().strip() != "Sig mismatch", "Error message is not what was expected"        

def test_file_sig_copy():
    """Test that streaming into a copy generates identical file"""
    test_fn = "test.txt"
    test_copy_fn = "test_copy.txt"
    with open(test_fn,"w") as out:
        out.write("FFUGJGJHGJGYNVNVvhgjgjyg76868uyu\n")
    test_sig = sig.file_sig(test_fn,out_copy=test_copy_fn)
    assert sig.file_sig_cmp(test_sig,inp_file=test_copy_fn), "Streamed file copy has a different signature from the original"
