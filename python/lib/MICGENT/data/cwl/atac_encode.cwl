class: CommandLineTool
cwlVersion: v1.0
hints:
  SoftwareRequirement:
    packages:
    - package: 'atac_encode'
      version:
      - '1.0'
baseCommand:
  - bds
inputs:
  - id: SampleID
    type: string
  - id: home
    type: string
  - id: config
    type: File
    inputBinding:
      prefix: '-c'
      position: 100
  - id: threads
    type: int?
    inputBinding:
      position: 0
      prefix: '-nth'
  - id: fastq1_1
    type: File
    inputBinding:
      position: 0
      prefix: '-fastq1_1'
  - id: fastq1_2
    type: File
    inputBinding:
      position: 0
      prefix: '-fastq1_2'
  - id: fastq2_1
    type: File?
    inputBinding:
      position: 0
      prefix: '-fastq2_1'
  - id: fastq2_2
    type: File?
    inputBinding:
      position: 0
      prefix: '-fastq2_2'
  - id: auto_detect_adapter
    type: boolean?
    inputBinding:
      position: 0
      prefix: '-auto_detect_adapter'
  - id: enable_idr
    type: boolean?
    inputBinding:
      position: 0
      prefix: '-enable_idr'
  - id: bds_home
    type: string
outputs:
  - id: out_json
    type: File
    outputBinding:
      glob: out/ENCODE_summary.json
  - id: out_archive
    type: File
    outputBinding:
      glob: $(inputs.SampleID + '.tar')
  - id: out_SampleID
    type: string
    outputBinding:
      outputEval: $(inputs.SampleID)

label: atac_encode
arguments:
  - position: 0
    valueFrom: $(inputs.home+"/atac.bds")
  - position: 200
    prefix: ''
    valueFrom: $('&& echo ' + inputs.SampleID + ' > out/SampleID && tar -cf ' + inputs.SampleID + '.tar out')
    separate: false
    shellQuote: false
requirements:
  - class: ShellCommandRequirement
  - class: InlineJavascriptRequirement
  - class: EnvVarRequirement
    envDef:
      PATH: $("\$PATH:"+inputs.bds_home)
      _JAVA_OPTIONS: '-Xms16G -Xmx16G -XX:ParallelGCThreads=1'
