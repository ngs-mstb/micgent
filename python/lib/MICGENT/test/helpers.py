"""Helper methods for testing.
This can be used with `import helpers` from test case files. This depends on the particular import rules
used by pytest, as described here: https://pytest.readthedocs.io/en/reorganize-docs/goodpractices.html
In particular, the import will most certainly break if __init__.py file is placed in the top `tests` directory.
"""
import pytest
import os
import re
from subprocess import call
import contextlib
import shutil

def get_conda_root(env=None):
    conda_root = None
    conda_exe = shutil.which("conda")
    if conda_exe:
        conda_root = os.path.dirname(os.path.dirname(conda_exe))
    return conda_root

skip_no_conda = pytest.mark.skipif(not get_conda_root(),
                    reason="requires active Conda environment")

def has_conda_env(conda_env):
    if get_conda_root():
        if call(["conda","activate",conda_env])==0 or call(["/bin/sh","-c","'source activate {}'".format(conda_env)])==0:
            return True
    return False

skip_no_conda_env_toil = pytest.mark.skipif(has_conda_env("toil"),
                    reason="requires Conda environment toil with Toil installation")


def grep(rx,fn):
    with open(fn,"rt",encoding="utf-8") as inp:
        s = inp.read()
        if re.search(rx,s,re.MULTILINE):
            return True
    return False

@contextlib.contextmanager
def mkchdir(dirname=None,create=True,exist_ok=True):
    curdir = os.getcwd()
    try:
        do_chdir = False
        if dirname is not None:
            if create:
                os.makedirs(dirname,exist_ok=exist_ok)
            dirname = abspath(dirname)
            do_chdir = True
        else:
            dirname = curdir
        if do_chdir:
            os.chdir(dirname)
        yield dirname
    finally:
        os.chdir(curdir)

def abspath(path):
    if not os.path.isabs(path):
        path = os.path.abspath(path)
    return path
