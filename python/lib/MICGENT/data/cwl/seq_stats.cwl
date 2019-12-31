#!/usr/bin/env cwl-runner

cwlVersion: v1.0

doc: Compute basic sequence statistics with Seqkit

requirements:
- class: InlineJavascriptRequirement

class: CommandLineTool

inputs:
  - id: inp_seqs
    type: File[]
    inputBinding:
      prefix:
      position: 100

  - id: out_basename
    type: string
    inputBinding:
      prefix: --out-file

  - id: threads
    type: int?
    inputBinding:
      prefix: --threads
      valueFrom: $(runtime.cores)

outputs:
  - id: out_stats
    type: File
    outputBinding:
      glob: $(out_basename)

baseCommand: ["seqkit", "--all"]
