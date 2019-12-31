"""Utility functions to massage a SAM file adopted from SRST2 package"""
from __future__ import unicode_literals
from __future__ import print_function
from __future__ import division
from __future__ import absolute_import
## Code adopted from SRST2 (https://github.com/katholt/srst2). Some of that
## code was contibuted earlier as patches to SRST2 by Andrey Tovchigrechko.
## The code was further modified here.
## The SRS2 license is BSD License.
## SRST2 Copyright is (c) 2013, Michael Inouye, Bernie Pope, Harriet Dashnow, Kathryn Holt.

from future import standard_library
standard_library.install_aliases()
from builtins import str
from builtins import *
from MGT.Logging import *
import re

def get_clips_cigar(cigar):
    """Return number and types of clipped bases on the left and right.
    @param cigar Cigar string (unparsed)
    @return (left_clip,right_clip), where each element is a dict
    with keys (type,length)"""
    #remove padding first if present
    #maybe padding is never present at the edges, but it is easier to just remove
    cigar = re.sub(r'\d+P','',cigar.strip())
    x = re.search(r'^(?P<length>\d+)(?P<type>[SH])',cigar)
    if x:
        left_clip = x.groupdict()
        left_clip["length"] = int(left_clip["length"])
    else:
        left_clip = dict(length=0,type=None)

    x = re.search(r'(?P<type>[SH])(?P<length>\d+)$',cigar)
    if x:
        right_clip = x.groupdict()
        right_clip["length"] = int(right_clip["length"])
    else:
        right_clip = dict(length=0,type=None)
    
    return (left_clip,right_clip)

def get_end_shift_cigar(cigar):
    """Return change in coordinate on the reference of the read end due to indels in CIGAR string"""
    shift = 0
    for edit_op in re.findall(r'(\d+)([ID])',cigar):
        shift += int(edit_op[0])*(-1 if edit_op[1] == 'I' else 1) 
    return shift

def get_unaligned_read_end_lengths_sam(fields,ref_len):
    """From SAM file line, compute clipped read length within reference.
    @param SAM fields as a sequence of strings
    @param ref_len length of reference for this SAM line
    @return tuple(left_clipped_len,right_clipped_len). If
    there are less than 10 fields in the record, return (0,0).
    """
    left_res = 0
    right_res = 0
    if len(fields) >= 10: 
        ## get (clipped) start position
        ali_clipped_start = int(fields[3])
        cigar = fields[5]
        ## get number and types of clipped bases on the left and right
        left_clip, right_clip = get_clips_cigar(cigar)
        left_res = min(ali_clipped_start,left_clip["length"])
        seq_start = ali_clipped_start
        if left_clip["type"] and left_clip["type"] == "S":
            seq_start -= left_clip["length"]
        ## get (hard-clipped) end position as start + len(seq)
        seq_hard_clipped_end = seq_start + len(fields[9]) + get_end_shift_cigar(cigar)
        ## seq end = hard end + right hard clip
        seq_end = seq_hard_clipped_end
        if right_clip["type"] and right_clip["type"] == "H":
            seq_end += right_clip["length"]
        ## aligned end = hard end - right soft clip
        ali_clipped_end = seq_hard_clipped_end
        if right_clip["type"] and right_clip["type"] == "S":
            ali_clipped_end -= right_clip["length"]
        ## right result = min(ref_len,right read end) - right aligned end
        right_res = min(ref_len,seq_end) - ali_clipped_end
    return (left_res,right_res)

def get_ref_length_sam(line,ref_lens):
    """Get reference length from @ LN tag and insert into dict"""
    if line.startswith('@SQ\t'):
        ref_search = re.search(r'\tSN:(\S+)\b',line)
        if ref_search:
            ref_name = ref_search.group(1)
            assert ref_name, "Empty reference name in {}".format(line)
            len_search = re.search(r'\tLN:(\d+)\b',line)
            assert len_search,"Could not find length tag in {}".format(line)
            ref_len = int(len_search.group(1))
            if ref_name in ref_lens:
                log.warning("Reference name is found second time in line {}".format(line))
            ref_lens[ref_name] = ref_len


def filter_bowtie_sam(raw_bowtie_sam,mod_bowtie_sam,max_mismatch=None,max_unaligned_overlap=10):
    """Fix sam flags for comprehensive pileup (multiple maytches per read) and filter out spurious alignments.
    The alignments where a read has too many mismatches or overhangs above the reference are removed.
    @param raw_bowtie_sam - SAM file as produced by Bowtie2 with options like 
    --very-sensitive-local --no-unal -a; in other words, potentially containing multiple
    matches to parts of the same read.
    @param max_mismatch maximum number of mismatches to tolerate [default is None - no filtering by that]
    @param max_unaligned_overlap maximum overhang of read over reference on either side to tolerate. The
    value to choose is application-specific - this filter will introduce uneven coverage around true
    insertions which are longer than this parameter value."""
    ref_lens = {}
    with open(raw_bowtie_sam,'r') as sam, open(mod_bowtie_sam, 'w') as sam_mod:
        for line in sam:
            if not line.startswith('@'):
                fields = line.split('\t')
                left_unali,right_unali = get_unaligned_read_end_lengths_sam(fields,ref_lens[fields[2].strip()])
                if left_unali > max_unaligned_overlap or right_unali > max_unaligned_overlap:
                    #log.debug("Excluding read from SAM file due to too long unaligned end overlapping the reference: {}".format(line))
                    continue
                flag = int(fields[1])
                flag = (flag - 256) if (flag & 256) else flag
                m = re.search("NM:i:(\d+)\s",line)
                if m != None:
                    num_mismatch = m.group(1)
                    if int(num_mismatch) <= int(max_mismatch):
                        sam_mod.write('\t'.join([fields[0], str(flag)] + fields[2:]))
                else:
                    log.info('Excluding read from SAM file due to missing NM (num mismatches) field: ' + fields[0])
                    num_mismatch = 0
            else:
                get_ref_length_sam(line,ref_lens)
                sam_mod.write(line)
            
    return ref_lens

## import package module and add argh entry points

def _main():
    import argh
    argh.dispatch_commands([
        filter_bowtie_sam
    ])


if __name__ == "__main__":
    _main()

