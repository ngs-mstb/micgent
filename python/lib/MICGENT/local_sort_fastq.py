"""Sort FASTQ query by reference when reodering is confined to a window"""

from . import local_sort
from . import util
from Bio.SeqIO.QualityIO import FastqGeneralIterator
import sys

def fastq_sort_by_ref_in_window(inp_ref=None,inp_que=None,
    window_size=1000,
    out=None,
    out_mode="w"):
    """Reorder inp_que iterator by inp_ref within window_size.

    This is the I/O interface to sort_by_ref_in_window.
    """
    assert not (inp_ref is None and inp_que is None),\
        "Reference and query cannot both come from the standard input stream"
    with util.as_stream_cm(inp_que) as inp_que,\
        util.as_stream_cm(inp_ref) as inp_ref,\
        util.as_stream_cm(out,out_mode) as out:
        inp_que = ((id,(seq,qual)) for (id, seq, qual) in 
                FastqGeneralIterator(inp_que))
        inp_ref = (_.rstrip() for _ in inp_ref)
        for id,(seq,qual) in local_sort.sort_by_ref_in_window(
            inp_ref=inp_ref,
            inp_que=inp_que,
            window_size=window_size):
            out.write(f"@{id}\n{seq}\n+\n{qual}\n")

## import package module and add argh entry points

def _main():
    import argh
    argh.dispatch_commands([
        fastq_sort_by_ref_in_window
    ])

if __name__ == "__main__":
    _main()
