from __future__ import absolute_import
from __future__ import unicode_literals
from __future__ import print_function
from __future__ import division
from future import standard_library
standard_library.install_aliases()
from builtins import str
from builtins import *

import os

def jinja_dump_template(env,name,out_fn,args=None):
    if not args:
        args = {}
    env.get_template(name).stream(**args).dump(out_fn)

def init_jinja(subdir):
    from jinja2 import Environment, PackageLoader, select_autoescape
    env = Environment(
        loader=PackageLoader('MICGENT', os.path.join('data/templates',subdir)),
        autoescape=select_autoescape(['html', 'xml'])
    )
    return env

