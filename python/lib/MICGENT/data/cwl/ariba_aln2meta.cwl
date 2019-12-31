#!/usr/bin/env cwl-runner
# This tool description was generated automatically by argparse2tool ver. 0.4.3-2
# To generate again: $ ariba --generate_cwl_tool
# Help: $ ariba --help_arg2cwl

cwlVersion: v1.0

class: CommandLineTool
baseCommand: ['ariba', '<command>', '<options>', 'aln2meta']

doc: |
  Converts multi-aln fasta and SNPs to metadata

inputs:
  
  genetic_code:
    type:
    - "null"
    - type: enum
      symbols: ['1', '4', '11']
    default: '11'
    doc: Number of genetic code to use. Currently supported 1,4,11 [%(default)s]
    inputBinding:
      prefix: --genetic_code 

  variant_only:
    type: ["null", boolean]
    default: False
    doc: Use this to flag all sequences as variant only. By default they are considered to be presence/absence
    inputBinding:
      prefix: --variant_only 

  aln_fasta:
    type: string
  
    doc: Multi-fasta file of alignments
    inputBinding:
      position: 1

  variants_tsv:
    type: string
  
    doc: TSV file of variants information
    inputBinding:
      position: 2

  coding_or_non:
    type:
      type: enum
      symbols: ['coding', 'noncoding']
    doc: Sequences are coding or noncoding. Must be one of - coding noncoding
    inputBinding:
      position: 3

  outprefix:
    type: string
  
    doc: Prefix of output filenames
    inputBinding:
      position: 4


outputs:
    []
