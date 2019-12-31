#!/usr/bin/env cwl-runner

cwlVersion: v1.0

class: CommandLineTool
baseCommand: ['python', '-m', 'MICGENT.gene_extractor', 'extract-contigs']
requirements:
  - class: InlineJavascriptRequirement
doc: |
  Extract contig sequences corresponding to Ariba report table
arguments:
  - valueFrom: $('report_out.tsv')
    position: 100
  - valueFrom: $('seq_out.fasta')
    position: 101
  - valueFrom: $('status_out.tsv')
    position: 102

inputs:
  one_seq_per_contig:
    type: boolean?
    default: false
    inputBinding:
      prefix: --one-seq-per-contig
      position: 1
  cut_to_ref:
    type: boolean?
    default: false
    inputBinding:
      prefix: --cut-to-ref
      position: 1      
  pad_assembled:
    type: int?
    default: 200
    inputBinding:
      prefix: --pad-assembled
      position: 1
  pad_gene:
    type: int?
    default: 200
    inputBinding:
      prefix: --pad-gene
      position: 1      
  sig_inp:
    type: string[]
    inputBinding:
      prefix: --sig-inp
      itemSeparator: ","
      position: 1      
  SampleID:
    type: string
    inputBinding:
      position: 10
  report:
    type: File
    inputBinding:
      position: 20
  assembled_genes:
    type: File
    inputBinding:
      position: 30
  assembled_seqs:
    type: File
    inputBinding:
      position: 40
  assemblies:
    type: File
    inputBinding:
      position: 50
  seq_map:
    type: File
    inputBinding:
      position: 60
  ariba_stdout:
    type: File
    inputBinding:
      position: 65
  ariba_stderr:
    type: File
    inputBinding:
      position: 70
  basecov_asm:
    type: File
    inputBinding:
      position: 75
  basecov_ref:
    type: File
    inputBinding:
      position: 80

outputs:
  report_out:
    type: File
    outputBinding:
      glob: report_out.tsv
  sequences:
    type: File
    outputBinding:
      glob: seq_out.fasta
  status_out:
    type: File
    outputBinding:
      glob: status_out.tsv
