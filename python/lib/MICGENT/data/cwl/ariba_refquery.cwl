#!/usr/bin/env cwl-runner
# This tool description was generated automatically by argparse2tool ver. 0.4.3-2
# To generate again: $ ariba --generate_cwl_tool
# Help: $ ariba --help_arg2cwl

cwlVersion: v1.0

class: CommandLineTool
baseCommand: ['ariba', 'refquery']

doc: |
  Get cluster or sequence info from prepareref output

inputs:
  
  prepareref_dir:
    type: string
  
    doc: Name of directory output by prepareref
    inputBinding:
      position: 1

  query_type:
    type:
      type: enum
      symbols: ['cluster', 'seq']
    doc: Use "cluster" to get the sequences in a cluster, or "seq" to get information about a sequence
    inputBinding:
      position: 2

  search_name:
    type: string
  
    doc: Name of cluster or sequence to search for
    inputBinding:
      position: 3


outputs:
    []
