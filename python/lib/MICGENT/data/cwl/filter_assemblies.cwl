#!/usr/bin/env cwl-runner

cwlVersion: v1.0

class: CommandLineTool
baseCommand: ['python', '-m', 'MICGENT.post_ariba', 'filter-assemblies']

doc: |
  Filter assemblies after the main extraction step

inputs:
  manifest:
    type: File
    inputBinding:
        prefix: --manifest
  sequences:
    type: File
    inputBinding:
        prefix: --contigs
  args:
    type: string?
    doc: Arguments for the filter as one-line YAML
    inputBinding:
      prefix: --args

outputs:
  manifest_out:
    type: File
    outputBinding:
      glob: manifest_out.tsv
  sequences_out:
    type: File
    outputBinding:
      glob: seq_out.fasta
  manifest_out_all:
    type: File
    outputBinding:
      glob: manifest_out_all.tsv
  manifest_out_sum:
    type: File
    outputBinding:
      glob: manifest_out_sum.tsv
  manifest_out_dict:
    type: File
    outputBinding:
      glob: manifest_out_dict.tsv
