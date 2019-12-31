#!/usr/bin/env cwl-runner
# This tool description was generated automatically by argparse2tool ver. 0.4.3-2
# To generate again: $ ariba --generate_cwl_tool
# Help: $ ariba --help_arg2cwl

cwlVersion: v1.0

class: CommandLineTool
baseCommand: ['ariba', 'pubmlstget']

doc: |
  Download species from PubMLST and make db

inputs:
  
  verbose:
    type: ["null", boolean]
    default: False
    doc: Be verbose
    inputBinding:
      prefix: --verbose 

  species:
    type: string
  
    doc: Species to download. Put it in quotes
    inputBinding:
      position: 1

  outdir:
    type: string
  
    doc: Name of output directory to be made (must not already exist)
    inputBinding:
      position: 2


outputs:
    []
