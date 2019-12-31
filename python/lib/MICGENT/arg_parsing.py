from __future__ import absolute_import
from __future__ import unicode_literals
from __future__ import print_function
from __future__ import division

from future import standard_library
standard_library.install_aliases()
from builtins import str
from builtins import *
from builtins import object

from . import util
from . import yaml_util
import sys
import argh
import collections
import functools
import argparse

def _join_sig_args(sig,args_passed):
    sig_args = sig.args if sig.args is not None else []
    sig_defaults = sig.defaults if sig.defaults is not None else []
    defs_dict = collections.OrderedDict(list(zip(sig_args[::-1], sig_defaults[::-1]))[::-1])
    assert len(sig_args) >= len(args_passed), \
        "Passed positional argument list is longer than method argument list"
    args_passed_dict = collections.OrderedDict(zip(sig_args,args_passed))
    return (args_passed_dict,defs_dict)

def update_sig_obj(wrapped,*args,**kwargs):
    #print("Wrapper received args:",args,kwargs)
    sig = argh.utils.get_arg_spec(wrapped)
    kwargs = kwargs.copy()
    ## take precedence over the config values, so do nothing
    ## for args that were passed positionally
    ##
    ## process defaults and keywords
    args_passed_dict, defs_dict = _join_sig_args(sig,args)
    for arg in sig.args:
        ## optional argument values passed as positional args always take precedence
        if arg not in args_passed_dict:
            #kwargs can override both positional and keyword remaining arguments
            #by setting them in the keyword dict
            if arg in kwargs:
                defs_dict[arg] = kwargs[arg]
    #print("Wrapper prepared args:",args_passed_dict,defs_dict)
    return dict(args_passed_dict=args_passed_dict,defs_dict=defs_dict)

def update_sig(wrapped,*args,**kwargs):
    o = update_sig_obj(wrapped=wrapped,*args,**kwargs)
    return functools.update_wrapper(functools.partial(
        wrapped,*list(o["args_passed_dict"].values()),**o["defs_dict"]),
        wrapped)

## Gotcha: parameter name 'function' trips something in
## a standard argh machinery, causes either exceptions or
## printing of a help message
def dump_sig(func,out_yaml=None,kwargs=None):
    """Save signature of the optional function arguments into YAML file (if provided).
    :param func: function or string with a full function import path
    :param kwargs: dict or YAML file name with keyword arguments to override defaults in func signature
    :return: OrderedDict (func name -> signature)
    """
    if util.is_string(func):
        func = util.import_name(func)
    if kwargs is None:
        kwargs = {}
    elif util.is_string(kwargs):
        kwargs = yaml_util.load_yaml(kwargs)
    o = update_sig_obj(wrapped=func,**kwargs)
    val = collections.OrderedDict([(func.__name__ , o["defs_dict"])])
    if out_yaml:
        yaml_util.dump_yaml(val,out_yaml)
    return val

class ConfigArgParser(argparse.ArgumentParser):
    def __init__(self,config_type=yaml_util.load_yaml, *args, **kwargs):
        super().__init__(add_help=False,
                         *args, **kwargs)
        self.add_argument("--config", default={}, type=config_type,
                          help="YAML config file with new default arguments for the command. \n"
                          "Should be a map of maps with command names as keys (dashes in all keys "
                          "must be converted to underscores). \n"
                          "Arguments from the command line will take precedence.")

class ArghParserChained(argh.ArghParser):
    """Subclass of argh.ArghParser that parses in two stages.
    The arguments loaded after the first stage can be used to construct
    further arguments for the parser API"""

    def __init__(self,first_parser,rest_arg="command",*args,**kwargs):
        self.rest_arg = rest_arg
        first_parser.add_argument("-h","--help",action="store_true",
                                  help="Print help and exit")
        first_parser.add_argument(rest_arg, nargs=argparse.REMAINDER,
                                  help="Command and its arguments")
        epilog = "This program can be called with additional arguments provided before the command as: \n"+\
            first_parser.format_help()
        super().__init__(epilog=epilog,*args, **kwargs)
        self.first_parser = first_parser

    def _parse_first(self,args=None):
        if args is None:
            args = sys.argv[1:]
        self.args = args
        first_args = self.first_parser.parse_args(args)
        self.rest_argv = getattr(first_args,self.rest_arg)
        delattr(first_args,self.rest_arg)
        self.first_args = first_args

    def dispatch(self, argv=None, *args, **kwargs):
        if argv is None:
            argv = self.rest_argv
        return super().dispatch(argv=argv,*args,**kwargs)

class ArghParserChainedConfig(ArghParserChained):
    """Subclass of ArghParserStaged that loads a config file in the first stage of parsing.
    The arguments loaded after the first stage can be used to construct
    further arguments for the parser API"""

    def __init__(self,first_parser=None, *args, **kwargs):
        if first_parser is None:
            first_parser = ConfigArgParser()
        super().__init__(first_parser,*args, **kwargs)

    def config_first(self):
        return self.first_args.config

    def add_commands(self, functions, argv=None, *args, **kwargs):
        self._parse_first(args=argv)
        config = self.config_first()
        ret = super().add_commands(
            [update_sig(f,**config.get(f.__name__,{})) for f in functions],
            *args,
            **kwargs)
        if self.first_args.help:
            self.print_help()
            self.exit(0)

## import package module and add argh entry points

def _main():
    #parser = ArghParserChainedConfig()
    parser = argh.ArghParser()

    parser.add_commands([
        dump_sig
    ])
    parser.dispatch()


if __name__ == "__main__":
    _main()
