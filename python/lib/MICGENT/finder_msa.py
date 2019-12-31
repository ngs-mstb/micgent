from __future__ import absolute_import
from __future__ import unicode_literals
from __future__ import print_function
from __future__ import division
from future import standard_library
standard_library.install_aliases()
from builtins import str
from builtins import *
from . import util

from subprocess import check_call

def add_seq_to_msa(
        seq_file,
        msa_file_inp,
        msa_file_out,
        flavor="addfragments",
        extra_mafft_args={}):
    args = util.get_arg_as_json(extra_mafft_args).copy()
    args[flavor] = seq_file
    cmd = ["mafft"]
    for key, val in list(args.items()):
        cmd += ["--{}".format(key),str(val)]
    cmd += [msa_file_inp]
    with open(msa_file_out,"w") as msa_out:
        check_call(cmd,stdout=msa_out)

def assign_to_ref_msa(
        msa,
        labels
):
    from . import seq_util
    mp = seq_util.ali_mapper(msa)
    import ipdb; ipdb.set_trace()

## import package module and add argh entry points

def _main():
    import argh
    argh.dispatch_commands([
        add_seq_to_msa,
        assign_to_ref_msa
    ])

if __name__ == "__main__":
    _main()
