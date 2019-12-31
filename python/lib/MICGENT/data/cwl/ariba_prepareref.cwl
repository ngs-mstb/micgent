#!/usr/bin/env cwl-runner
# This tool description was generated automatically by argparse2tool ver. 0.4.3-2
# To generate again: $ ariba --generate_cwl_tool
# Help: $ ariba --help_arg2cwl

cwlVersion: v1.0

class: CommandLineTool
baseCommand: ['ariba', 'prepareref']

doc: |
  Prepare reference data for input to "run". REQUIRED: -f/--fasta, and also either -m/--metadata or --all_coding must be used

inputs:
  
  fasta_files:
    type: 
      type: array
      items: File
      inputBinding:
        prefix: --fasta 
    doc: REQUIRED. Name of fasta file. Can be used more than once if your sequences are spread over more than on file

  tsv_files:
    type:
    - "null" 
    - type: array
      items: File
      inputBinding:
        prefix: --metadata 
    doc: Name of tsv file of metadata about the input sequences. Can be used more than once if your metadata is spread over more than one file. Incompatible with --all_coding

  all_coding:
    type:
    - "null"
    - type: enum
      symbols: ['yes', 'no']
    doc: Use this if you only have a fasta of presence absence sequences as input, and no metadata. Use "yes" if all sequences are coding, or "no" if they are all non-coding. Incompatible with -m/--metadata
    inputBinding:
      prefix: --all_coding 

  no_cdhit:
    type: ["null", boolean]
    default: False
    doc: Do not run cd-hit. Each input sequence is put into its own "cluster". Incompatible with --cdhit_clusters.
    inputBinding:
      prefix: --no_cdhit 

  cdhit_clusters:
    type: ["null", File]
    doc: File specifying how the sequences should be clustered. Will be used instead of running cdhit. Format is one cluster per line. Sequence names separated by whitespace. Incompatible with --no_cdhit
    inputBinding:
      prefix: --cdhit_clusters 

  cdhit_min_id:
    type: ["null", float]
    default: 0.9
    doc: Sequence identity threshold (cd-hit option -c) [%(default)s]
    inputBinding:
      prefix: --cdhit_min_id 

  cdhit_min_length:
    type: ["null", float]
    default: 0.0
    doc: length difference cutoff (cd-hit option -s) [%(default)s]
    inputBinding:
      prefix: --cdhit_min_length 

  min_gene_length:
    type: ["null", int]
    default: 6
    doc: Minimum allowed length in nucleotides of reference genes [%(default)s]
    inputBinding:
      prefix: --min_gene_length 

  max_gene_length:
    type: ["null", int]
    default: 10000
    doc: Maximum allowed length in nucleotides of reference genes [%(default)s]
    inputBinding:
      prefix: --max_gene_length 

  genetic_code:
    type:
    - "null"
    - type: enum
      symbols: ["1", "4", "11"]
    default: "11"
    doc: Number of genetic code to use. Currently supported 1,4,11 [%(default)s]
    inputBinding:
      prefix: --genetic_code 

  force:
    type: ["null", boolean]
    default: False
    doc: Overwrite output directory, if it already exists
    inputBinding:
      prefix: --force 

  threads:
    type: ["null", int]
    default: 1
    doc: Number of threads (currently only applies to cdhit) [%(default)s]
    inputBinding:
      prefix: --threads 

  verbose:
    type: ["null", boolean]
    default: False
    doc: Be verbose
    inputBinding:
      prefix: --verbose 

  outdir:
    type: string
  
    doc: Output directory (must not already exist)
    inputBinding:
      position: 1


outputs:
    []
