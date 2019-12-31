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
from . import arg_parsing
import argh
import hmac
import sys
import os


class SignatureError(Exception):
    pass

## For speed, we use native
## keyed hashing in blake2b instead of through hmac:
## https://docs.python.org/3/library/hashlib.html#hashlib.blake2b
## See also this post for comparison with different flavors of SHA:
## https://research.kudelskisecurity.com/2017/03/06/why-replace-sha-1-with-blake2/
## Blake is only available starting in Python 3.6
def str_sig(msg, key="123", codec="utf8", digestmod="blake2b"):
    """Compute hex digest of a file mixing with the (secret) key.
    To keep key really secret when calling with CLI, use the config file interface for passing the key.

    Accepts either str or byte/bytearray objects.
    """
    if util.is_string(key):
        key = key.encode(codec)
    if util.is_string(msg):
        msg = msg.encode(codec)
    if digestmod == "blake2b":
        from hashlib import blake2b
        dig = blake2b(key=key)
        dig.update(msg)
        ret = dig.hexdigest()
    else:
        ret = hmac.new(key, msg=msg, digestmod=digestmod).hexdigest()
    return ret


def file_sig(inp_file, key="123", chunk_size=2 ** 20, codec="utf8", digestmod="blake2b", out_file=None, extra_data="", out_copy=None):
    """Compute hex digest of a file mixing with the (secret) key.
    To keep key really secret when calling with CLI, use the config file interface for passing the key.
    Hint: use YAML anchors to propagate the same key to function names inside the same YAML config

    :param out_file: If defined, write signature to this file name instead of returning the value. Containing directories
    will be created if not already exist (for ease of integration with Galaxy).
    :param extra_data: Optional extra string or byte array that will be used as the first block to hash. It plays the
    same role as the `person` argument to blake2 algorithm, but does not have the latter length restriction of 16 bytes.
    :param out_copy: If defined, stream a copy of the input file data into this file name, as the input file is being 
    read for signature generation. This addresses a use case where a signature of a large file has to be generated
    while creating a protected copy of the input file in a single transaction and avoiding reading the file data twice.

    For example, if the workflow records a signature of a FASTQ file from a NFS before computing anything on the file,
    it makes sense to create a local copy that will then be used for computations. Otherwise, it is theoretically possible
    that the input file on the NFS will be modified externally between the moment the signature is created and the
    moment the file is used in computations.
    """
    if util.is_string(key):
        key = key.encode(codec)
    if util.is_string(extra_data):
        extra_data = extra_data.encode(codec)        
    if digestmod == "blake2b":
        from hashlib import blake2b
        hm = blake2b(key=key)
    else:    
        hm = hmac.new(key, digestmod=digestmod)
    ## the blake2 implementation seems to be making no difference between key, salt and person arguments -
    ## it treats them together simply as a data block for the initial update:
    ## https://github.com/dchest/pyblake2/search?q=personal&unscoped_q=personal
    ## It XORs this data block with the 32 byte hash state
    ## https://github.com/dchest/pyblake2/blob/master/impl/blake2s.c
    ## and it does it only once at the beginning. blake2s_update() in blake2s.c does not access the parameter structure anymore.
    ## Therefore, it seems reasonable to use extra extra_data string as the first data block, instead of trying to put it into
    ## salt or person arguments, which are limited to 16 bytes each. It is a different question though if having extra_data longer than
    ## 16 bytes affects the strength of the digest.
    if extra_data:
        hm.update(extra_data)
    buf = bytearray(chunk_size)
    try:
        if out_copy:
            out_copy = open(out_copy,"wb")
        with open(inp_file, "rb") as inp:
            while True:
                n_read = inp.readinto(buf)
                if n_read:
                    buf_read = buf[:n_read]
                    hm.update(buf_read)
                    if out_copy:
                        out_copy.write(buf_read)
                else:
                    break
    finally:
        if out_copy:
            out_copy.close()
    dig = hm.hexdigest()
    if not out_file:
        return dig
    else:
        util.make_pardir(out_file)
        with open(out_file,"w") as out:
            out.write(dig)

def file_sig_cmp(sig, inp_file, key="123", chunk_size=2 ** 20, codec="utf8", digestmod="blake2b", extra_data=""):
    """Compute hex digest of a file mixing with a secret key, and compare with a previously compared signature.
    To keep key really secret when calling with CLI, use the config file interface for passing the key.
    Hint: use YAML anchors to propagate the same key to function names inside the same YAML config

    :param extra_data: Optional extra string or byte array that will be used as the first block to hash. It plays the
    same role as the `person` argument to blake2 algorithm, but does not have the latter length restriction of 16 bytes.

    :return: True or False. Note that in argh CLI mode, it will print True or False string, but the exit code will
    still be 0.
    """
    sig_f = file_sig(inp_file=inp_file, key=key, chunk_size=chunk_size, codec=codec, digestmod=digestmod, extra_data=extra_data)
    return hmac.compare_digest(sig_f, sig)

## Comment wrap_errors out until https://github.com/neithere/argh/issues/118 is fixed
#@argh.wrap_errors([SignatureError],processor=lambda excinfo: '{0}'.format(excinfo))
def file_sig_cmp_msg(sig_file, inp_file, msg="File signature mismatch", key="123", chunk_size=2 ** 20, codec="utf8", digestmod="blake2b", extra_data=""):
    """Version of file_sig_cmp intended for use in CLI.
    It exits with code 1 and a message in stderr"""
    args = locals()
    del args["msg"]
    del args["sig_file"]
    if not os.path.isfile(sig_file):
        sys.exit("Required signature file does not exist")
    max_sig_size = 2 ** 20
    try:
        with open(sig_file,"r") as inp:
            args["sig"] = inp.read(max_sig_size).strip()
    except:
        sys.exit("Unable to read signature file")
    if not file_sig_cmp(**args):
        sys.exit(msg)

def bioseq_sig_iter(inp_file, key="123", codec="utf8", digestmod="blake2b", seq_format="fasta", with_file_name=False):
    from Bio import SeqIO
    inp = None
    try:
        if inp_file.endswith(".gz"):
            inp = util.open_gzip_text23(inp_file, "rt")
        else:
            inp = open(inp_file, "rt")
        for record in SeqIO.parse(inp, seq_format):
            out_rec = (record.id, str_sig(str(record.seq), key=key, codec=codec, digestmod=digestmod))
            if with_file_name:
                out_rec = (inp_file,) + out_rec
            yield out_rec
    finally:
        if inp is not None:
            inp.close()

def bioseq_sig(inp_file, key="123", codec="utf8", digestmod="blake2b", seq_format="fasta", with_file_name=False,
               out_csv="-", out_sep="\t", out_codec="utf8"):
    from Bio import SeqIO
    out = None
    try:
        if out_csv == "-":
            out = sys.stdout
        else:
            out = open(out_csv,mode="wt",encoding=out_codec)
        for out_rec in bioseq_sig_iter(inp_file=inp_file,
            key=key,
            codec=codec,
            digestmod=digestmod,
            seq_format=seq_format,
            with_file_name=with_file_name):
            out.write(out_sep.join(out_rec) + "\n")
    finally:
        if out is not None and out_csv != "-":
            out.close()


def dir_sig(inp_dir, key="123",
            chunk_size=2 ** 20,
            codec="utf8",
            digestmod="blake2b",
            with_dir_name=False,
            out_csv=None,
            out_sep="\t",
            out_codec="utf8",
            no_recurse=False,
            filt=None,
            filt_re_format="glob",
            filt_re_flags=0,
            filt_re_apply="basename"
            ):
    """Walk a directory tree and compute hex digests for all files
    """
    res = []
    out = None
    try:
        if out_csv:
            if out_csv == "-":
                out = sys.stdout
            else:
                out = open(out_csv, mode="wt", encoding=out_codec)
        for inp_start,inp_end in util.glob_files(files_globs=inp_dir,
                                        filt=filt,
                                        filt_re_format=filt_re_format,
                                        filt_re_flags=filt_re_flags,
                                        filt_re_apply=filt_re_apply,
                                        recurse=not no_recurse,
                                        return_split_path=True,
                                        entry_type="f"
                                        ):
            inp_file = os.path.join(inp_start,inp_end)
            if not os.path.islink(inp_file):
                out_rec = (inp_end,
                           file_sig(inp_file=inp_file, key=key,
                                    chunk_size=chunk_size, codec=codec, digestmod=digestmod))
                if with_dir_name:
                    out_rec = (inp_start,) + out_rec
                if not out:
                    res.append(out_rec)
                else:
                    out.write(out_sep.join(out_rec) + "\n")
    finally:
        if out is not None and out_csv != "-":
            out.close()
    if not out_csv:
        return res


## import package module and add argh entry points

def _main():
    from . import arg_parsing
    parser = arg_parsing.ArghParserChainedConfig()
    parser.add_commands([
        str_sig,
        file_sig,
        file_sig_cmp,
        file_sig_cmp_msg,
        bioseq_sig,
        dir_sig
    ])
    parser.dispatch()


if __name__ == "__main__":
    _main()
