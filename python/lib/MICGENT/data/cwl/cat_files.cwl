class: CommandLineTool
cwlVersion: v1.0
doc: >-
  Concatenate several files together. This currently uses Unix `cat` utility,
  and should work correctly on gzip files.
baseCommand:
  - cat
inputs:
  - id: files
    type:
      type: array
      items: File
      inputBinding:
        prefix: ''
        separate: false
    inputBinding:
      position: 1
  - id: out_nameroot
    type: string?
    default: 'out'
outputs:
  - id: output
    type: File
    outputBinding:
      glob: '$(inputs.out_nameroot+inputs.files[0].nameext)'
arguments:
  - position: 100
    prefix: ''
    shellQuote: false
    valueFrom: '$("> "+inputs.out_nameroot+inputs.files[0].nameext)'
requirements:
  - class: InlineJavascriptRequirement
  - class: ShellCommandRequirement
