#!/usr/bin/env cwl-runner
# This tool description was generated automatically by argparse2tool ver. 0.4.3-2
# To generate again: $ ariba --generate_cwl_tool
# Help: $ ariba --help_arg2cwl

cwlVersion: v1.0

class: CommandLineTool
baseCommand: ['ariba', 'expandflag']

doc: |
  Expands flag column of report file

inputs:
  
  infile:
    type: string
  
    doc: Name of input report TSV file
    inputBinding:
      position: 1

  outfile:
    type: string
  
    doc: Name of output report TSV file
    inputBinding:
      position: 2


outputs:
    []
