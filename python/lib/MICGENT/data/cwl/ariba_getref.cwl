#!/usr/bin/env cwl-runner
# This tool description was generated automatically by argparse2tool ver. 0.4.3-2
# To generate again: $ ariba --generate_cwl_tool
# Help: $ ariba --help_arg2cwl

cwlVersion: v1.0

class: CommandLineTool
baseCommand: ['ariba', 'getref']

doc: |
  Download reference data

inputs:
  
  debug:
    type: ["null", boolean]
    default: False
    doc: Do not delete temporary downloaded files
    inputBinding:
      prefix: --debug 

  version:
    type: ["null", string]
    doc: Version of reference data to download. If not used, gets the latest version. Only applies to card and megares
    inputBinding:
      prefix: --version 

  db:
    type:
      type: enum
      symbols: ['argannot', 'card', 'megares', 'plasmidfinder', 'resfinder', 'srst2_argannot', 'vfdb_core', 'vfdb_full', 'virulencefinder']
    doc: Database to download. Must be one of - argannot card megares plasmidfinder resfinder srst2_argannot vfdb_core vfdb_full virulencefinder
    inputBinding:
      position: 1

  outprefix:
    type: string
  
    doc: Prefix of output filenames
    inputBinding:
      position: 2


outputs:
    []
