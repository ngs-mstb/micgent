#!/usr/bin/env cwl-runner
# This tool description was generated automatically by argparse2tool ver. 0.4.3-2
# To generate again: $ ariba --generate_cwl_tool
# Help: $ ariba --help_arg2cwl

cwlVersion: v1.0

class: CommandLineTool
baseCommand: ['ariba', 'summary']

doc: |
  Summarise multiple reports made by "run". Files must be listed after the output file and/or the option --fofn must be used. If both used, all files in the filename specified by --fofn AND the files listed after the output file will be used as input.

inputs:
  
  fofn:
    type: File?
    doc: File of filenames of ariba reports to be summarised. Must be used if no input files listed after the outfile. The first column should be the filename. An optional second column can be used to specify a sample name for that file, which will be used instead of the filename in output files. Columns separated by whitespace.
    inputBinding:
      prefix: --fofn 

  preset:
    type:
    - "null"
    - type: enum
      symbols: ['minimal', 'cluster_small', 'cluster_all', 'cluster_var_groups', 'all', 'all_no_filter']
    doc: Shorthand for setting --cluster_cols,--col_filter,--row_filter,--v_groups,--variants. Using this overrides those options
    inputBinding:
      prefix: --preset 

  cluster_cols:
    type: ["null", string]
    default: match
    doc: Comma separated list of cluster columns to include. Choose from - assembled, match, ref_seq, pct_id, ctg_cov, known_var, novel_var [%(default)s]
    inputBinding:
      prefix: --cluster_cols 

  col_filter:
    type:
    - "null"
    - type: enum
      symbols: ['y', 'n']
    default: y
    doc: Choose whether columns where all values are "no" or "NA" are removed [%(default)s]
    inputBinding:
      prefix: --col_filter 

  no_tree:
    type: ["null", boolean]
    default: False
    doc: Do not make phandango tree
    inputBinding:
      prefix: --no_tree 

  row_filter:
    type:
    - "null"
    - type: enum
      symbols: ['y', 'n']
    default: y
    doc: Choose whether rows where all values are "no" or "NA" are removed [%(default)s]
    inputBinding:
      prefix: --row_filter 

  min_id:
    type: ["null", float]
    default: 90
    doc: Minimum percent identity cutoff to count as assembled [%(default)s]
    inputBinding:
      prefix: --min_id 

  only_clusters:
    type: ["null", string]
    doc: Only report data for the given comma-separated list of cluster names, eg - cluster1,cluster2,cluster42
    inputBinding:
      prefix: --only_clusters 

  v_groups:
    type: ["null", boolean]
    default: False
    doc: Show a group column for each group of variants
    inputBinding:
      prefix: --v_groups 

  known_variants:
    type: ["null", boolean]
    default: False
    doc: Report all known variants
    inputBinding:
      prefix: --known_variants 

  novel_variants:
    type: ["null", boolean]
    default: False
    doc: Report all novel variants
    inputBinding:
      prefix: --novel_variants 

  verbose:
    type: ["null", boolean]
    default: False
    doc: Be verbose
    inputBinding:
      prefix: --verbose 

  outprefix:
    type: string
  
    doc: Prefix of output files
    inputBinding:
      position: 1

  infiles:
    type:
    - "null"
    - type: array
      items: File
  
    doc: Files to be summarised
    inputBinding:
      position: 2


outputs:
    []
