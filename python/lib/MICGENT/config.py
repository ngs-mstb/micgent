### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##
#
#   See COPYING file distributed along with the MDDV package for the
#   copyright and license terms.
#
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##

"""Methods for loading config files"""
from __future__ import absolute_import
from __future__ import unicode_literals
from __future__ import print_function
from __future__ import division

from future import standard_library
standard_library.install_aliases()
from builtins import str
from builtins import *
from builtins import object
from . import resources
from . import util
import os


def env_vars_to_conf(env,prefix,low_keys=True):
    """Convert env vars like PREFIX_KEY=value to key=value"""
    out = {}
    prefix = prefix + "_"
    for key, val in list(env.items()):
        if key.startswith(prefix):
            key = key[len(prefix):]
            if len(key):
                out[key.lower()] = val
    return out

def load_config(config_file=None,
                vars={},
                vars_default={},
                pkg=None,
                prefix_env=None,
                config_file_base=None,
                config_var_sfx="CONFIG"
                ):
    """Try to load jsonnet config file from a variety of default locations.
    The order of preference: config_file; location set by an environment variable;
    current directory; package data directory.
    The base file name to look for in current directory and package data is given by
    config_file_base. If None, then it is constructed from prefix_env in lower case.
    If prefix_env is None, it is constructed from the last component of pkg in upper case.
    In other words, the following locations will be checked if a call is
    load_config(pkg="My.Pkg",vars=vars):

        file pointed to by MY_PKG_CONFIG (note dot replaced with _; CONFIG suffix can be
        overridden by config_var_sfx)
        mypkg.jsonnet in the current dir
        mypkg.jsonnet installed inside Python package mypkg at directory data

    Once config file is identified, the variables for substitution will be prepared:
    All environment variables that look like MYPKG_XXX=yyy are converted to xxx=yyy and
    updated with dictionary vars; the resulting dictionary is used to update
    vars_default; vars_default is passed to jsonnet for variable substitution wherever
    std.extVar("xxx") is used in the config file. Thus, values in vars override any other
    values, while values in vars_default are used as defaults before anything else is set.
    Return value if dictionary produced by jsonnet.evaluate_snippet().
    """
    import _jsonnet
    import json

    if config_file is None and pkg is None and prefix_env is None:
        raise ValueError("Either config_file or pkg or prefix_env argument must be defined")

    if prefix_env is None:
        if pkg is not None:
            prefix_env = pkg.replace(".","_").upper()

    config_var_sfx = config_var_sfx.lower()

    all_vars = vars_default.copy()

    if prefix_env is not None:
        all_vars.update(env_vars_to_conf(os.environ,prefix=prefix_env))

    all_vars.update(vars)

    if config_file is None:
        config_file = all_vars.get(config_var_sfx,None)

    if config_file is None:
        if config_file_base is None:
            config_file_base = prefix_env.lower() + ".jsonnet"
        config_file = config_file_base

        ## check the local dir
        if os.path.exists(config_file):
            with util.open_text_py23(config_file,"r") as _:
                config_str = _.read()
        elif pkg is not None:
            config_str = resources.get_pkg_data_string(pkg=pkg,name=config_file)
    else:
        with util.open_text_py23(config_file,"r") as _:
            config_str = _.read()


    conf = json.loads(str(_jsonnet.evaluate_snippet(config_file, config_str, ext_vars=all_vars)))

    return conf

def save_config_json(config,config_file):
    util.save_json(config,config_file)

def _jsonnet_vars_convert(o):
    return dict( (key,",".join(str(_) for _ in val) if util.is_sequence(val) else str(val) ) for (key,val) in list(o.items()))

"""Class to load JSONNET config once, and interpret many times"""
class config_interpreter(object):

    def __init__(self,config_file, vars_default={}):
        with util.open_text_py23(config_file,"r") as _:
            self.config_str = _.read()
        self.config_file = config_file
        self.vars_default = vars_default

    def __call__(self,vars={}):
        import _jsonnet
        import json

        all_vars = self.vars_default.copy()
        all_vars.update(vars)
        all_vars = _jsonnet_vars_convert(all_vars)
        # self.config_file here is needed only to annotate stacktraces in case of errors
        return json.loads(_jsonnet.evaluate_snippet(self.config_file, self.config_str, ext_vars=all_vars))

def load_subconfig_interpreter(obj,key_config="config",key_vars="vars",vars_default={}):
    vars_def=obj.get(key_vars,{}).copy()
    vars_def.update(vars_default)
    return config_interpreter(config_file=obj[key_config],vars_default=vars_def)

