#!/usr/bin/env cwl-runner

cwlVersion: v1.0

class: CommandLineTool
baseCommand: ['python', '-m', 'MICGENT.sig']
hints:
  - class: SoftwareRequirement
    packages:
    - package: 'ngs-mstb'
      version:
      - '1.0'
requirements:
  - class: InlineJavascriptRequirement
doc: |
  Compute cryptographic signature of a file and optionally create a copy of the file streamed while computing the signature

arguments:
  - valueFrom: $('file-sig')
    position: 2
  - valueFrom: $('--out-file')
    position: 3
  - valueFrom: $('sig.txt')
    position: 4

inputs:
  config:
    type: File?
    inputBinding:
      prefix: --config 
      position: 1
  out_copy_prefix:
    type: string?
    inputBinding:
      position: 10
      prefix: "--out-copy"
      valueFrom: >
        ${ if(self) 
          { return self+inputs.inp_file.basename; }
          else
          { return null; }
        }
  key:
    type: string?
    inputBinding:
      position: 10
      prefix: "--key"        
  inp_file:
    type: File
    inputBinding:
      position: 20

outputs:
  sig:
    type: string
    outputBinding:
      glob: sig.txt
      loadContents: true
      outputEval: $(self[0].contents)    
  out_copy:
    type: File?
    outputBinding:
      glob: $(inputs.out_copy_prefix+inputs.inp_file.basename)
