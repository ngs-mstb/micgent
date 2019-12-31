#!/usr/bin/env cwl-runner

class: CommandLineTool
cwlVersion: v1.0

doc: Use BBnorm from BBMap toolset for digital normalization of read depth
hints:
- class: ResourceRequirement
  ramMin: 64000
requirements:
- class: InlineJavascriptRequirement
  expressionLib:
    $import: micgent_js_lib.yaml
- class: SchemaDefRequirement
  types:
    $import: bbduk_types.yaml
- class: EnvVarRequirement
  envDef:
    #still need this because if comes as somehow defined (it should not), it will overwrite
    #whatever is passed to the BBTools wrapper script
    _JAVA_OPTIONS: "-Xms$(runtime.ram)m -Xmx$(runtime.ram)m -XX:ParallelGCThreads=1"

inputs:
  - id: SampleID
    type: string

  - id: prefix?
    type: string
    default: "_nr"

  - id: target
    type: int?
    inputBinding:
      position: 8
      prefix: "target="
      separate: false
    default: 200

  - id: mindepth
    type: int?
    inputBinding:
      position: 8
      prefix: "mindepth="
      separate: false
    default: 2

  - id: passes
    type: int?
    inputBinding:
      position: 8
      prefix: "passes="
      separate: false

  - id: fixspikes
    type:  boolean?
    inputBinding:
      position: 10
      prefix: "fixspikes"

  - id: ecc
    type:  boolean?
    inputBinding:
      position: 10
      prefix: "ecc"

  - id: k
    type: int?
    inputBinding:
      position: 12
      prefix: "k="
      separate: false

  - id: threads
    type: int?
    inputBinding:
      position: 16
      prefix: "threads="
      separate: false
      valueFrom: $(runtime.cores)

  - id: inp_seq1
    type: File
    inputBinding:
      position: 18
      prefix: "in="
      separate: false

  - id: inp_seq2
    type: File?
    inputBinding:
      position: 20
      prefix: "in2="
      separate: false

outputs:

  - id: out_seqs
    type: File[]
    outputBinding:
      glob: '$(inputs.SampleID)$(inputs.prefix)_[0-9].fastq.gz'
      outputEval: $(sortFileObjectsByName(self))

  - id: out_SampleID
    type: string
    outputBinding:
      outputEval: $(inputs.SampleID)

  - id: out_stats
    type: stdout

arguments:

  - valueFrom: $("out="+inputs.SampleID+inputs.prefix+"_1.fastq.gz")
  - valueFrom: $("out2="+inputs.SampleID+inputs.prefix+"_2.fastq.gz")
  - valueFrom: $("-Xmx"+runtime.ram+"m")

baseCommand: ["bbnorm.sh"]
stdout: $(inputs.SampleID+inputs.prefix+".bbnorm.stats")