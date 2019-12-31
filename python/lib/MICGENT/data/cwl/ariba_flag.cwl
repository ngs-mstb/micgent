#!/usr/bin/env cwl-runner
# This tool description was generated automatically by argparse2tool ver. 0.4.3-2
# To generate again: $ ariba --generate_cwl_tool
# Help: $ ariba --help_arg2cwl

cwlVersion: v1.0

class: CommandLineTool
baseCommand: ['ariba', 'flag']

doc: |
  Translate the meaning of a flag

inputs:
  
  flag_in:
    type: int
  
    doc: Flag to be translated (an integer)
    inputBinding:
      position: 1


outputs:
    []
