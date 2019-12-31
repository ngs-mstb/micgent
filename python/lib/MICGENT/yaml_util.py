"""Helper function for YAML
"""
from __future__ import absolute_import
from __future__ import unicode_literals
from __future__ import print_function
from __future__ import division
from future import standard_library
standard_library.install_aliases()
from builtins import str
from builtins import *

from . import util
from ruamel import yaml
import codecs
import collections

## copied from scriptcwl
def is_multiline(s):
    """Return True if a str consists of multiple lines.

    Args:
        s (str): the string to check.

    Returns:
        bool
    """
    return len(s.splitlines()) > 1

## copied from scriptcwl
def str_presenter(dmpr, data):
    """Return correct str_presenter to write multiple lines to a yaml field.

    Source: http://stackoverflow.com/a/33300001
    """
    if is_multiline(data):
        return dmpr.represent_scalar('tag:yaml.org,2002:str', data, style='|')
    return dmpr.represent_scalar('tag:yaml.org,2002:str', data)

def ordered_dict_presenter(dmpr, data):
    """Suppress saving OrderedDict as !!omap, save a regular map instead
    """
    ## this function already receives dict.items() when called from RoundTripRepresenter
    return dmpr.represent_mapping('tag:yaml.org,2002:map', data)

## The lines below modify the behaviour of ruamel globally

yaml.add_representer(str, str_presenter)
yaml.add_representer(collections.OrderedDict, ordered_dict_presenter)

yaml.representer.RoundTripRepresenter.add_representer(str, str_presenter)
yaml.representer.RoundTripRepresenter.add_representer(collections.OrderedDict, ordered_dict_presenter)

def dumps_yaml(x):
    return yaml.dump(x, Dumper=yaml.RoundTripDumper,default_flow_style=False)

def dump_yaml(x,fname,encoding = 'utf-8',header=None):

    with codecs.open(fname, 'wb', encoding=encoding) as yaml_file:
        if header:
            yaml_file.write('{}\n'.format(header))
        yaml_file.write(dumps_yaml(x))

def loads_yaml(inp):
    return yaml.load(inp, yaml.RoundTripLoader)

def load_yaml(fname,encoding = 'utf-8'):
    with codecs.open(fname, 'rb', encoding=encoding) as yaml_file:
        return loads_yaml(yaml_file)

def to_yaml(x):
    """Take YAML-compatible data structure, and return Ruamel data structure from round-trip loader.

    The return value then can be used to add YAML comments or call other YAML-specific Ruamel object calls.
    """
    return loads_yaml(dumps_yaml(x))

def get_arg_as_yaml(x):
    if util.is_string(x):
        x = yaml.safe_load(x)
    return x

def diff_yaml(x,y):
    """Recursively compare two YAML-compatible data structures and return parts that are different.

    This is utility method to use in testing results of transforming such deeply recursive structures
    against expected results.
    """
    ret = None
    if util.is_sequence(x):
        if util.is_sequence(y):
            ret = []
            for x1,y1 in zip(x,y):
                if x1!=y1: 
                    ret.append(diff_yaml(x1,y1))
    elif util.is_mapping(x):
        if util.is_mapping(y):
            ret = {}
            if x != y:
                for xkey,xval in x.items():
                    if not xkey in y:
                        ret[xkey] = collections.OrderedDict(Comparison=False,x=xval,y=None)
                    else:
                        yval = y[xkey]
                        if xval != yval:
                            ret[xkey] = diff_yaml(xval,yval)
                for ykey,yval in y.items():
                    if not ykey in x:
                        ret[ykey] = collections.OrderedDict(Comparison=False,x=None,y=yval)
    if ret is None:
        ret = collections.OrderedDict(Comparison=(x==y),x=x,y=y)
    return ret


__all__ = ['yaml','dumps_yaml','dump_yaml','loads_yaml','load_yaml','get_arg_as_yaml','is_multiline','diff_yaml']
