#!/usr/bin/env cwl-runner

class: CommandLineTool
cwlVersion: v1.0

doc: Use BBDuk from BBMap toolset to trim adapters and primers

hints:
- class: ResourceRequirement
  ramMin: 4096
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
  # default k, mink and hdist are taken from bbduk author's post on adapter trimming
  # https://www.biostars.org/p/155165/#268947
  # Except that we set mink=2 in order to catch partial read-through adapters
  # at the expense of losing some small amount of legitimate end fragments (which
  # are usually low quality anyway). hdist2 is set to 0.
  - id: SampleID
    type: string

  - id: prefix?
    type: string
    default: ""

  - id: tpe
    type: boolean?
    inputBinding:
      position: 10
      prefix: "tpe"

  - id: tbo
    type: boolean?
    inputBinding:
      position: 10
      prefix: "tbo"
    default: true

  - id: copyundefined
    type: boolean?
    inputBinding:
      position: 10
      prefix: "copyundefined"
    default: true

  - id: mm
    type:  bbduk_types.yaml#flag?
    inputBinding:
      position: 10
      prefix: "mm="
      separate: false

  - id: k
    type: int?
    inputBinding:
      position: 12
      prefix: "k="
      separate: false
    default: 23

  - id: mink
    type: int?
    inputBinding:
      position: 12
      prefix: "mink="
      separate: false
    default: 2

  - id: hdist
    type: int?
    inputBinding:
      position: 14
      prefix: "hdist="
      separate: false
    default: 1

  - id: hdist2
    type: int?
    inputBinding:
      position: 14
      prefix: "hdist2="
      separate: false
    default: 0

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

  - id: adapter_literals
    type: string[]?
    inputBinding:
      position: 22
      prefix: "literal="
      itemSeparator: ","
      separate: false

  - id: adapter_file
    type: File?
    inputBinding:
      position: 24
      prefix: "ref="
      separate: false

  - id: ktrim
    type: bbduk_types.yaml#ktrim
    inputBinding:
      position: 26
      prefix: "ktrim="
      separate: false

  - id: rcomp
    type:  bbduk_types.yaml#flag?
    inputBinding:
      position: 28
      prefix: "rcomp="
      separate: false

  - id: samplerate
    type: float?
    inputBinding:
      position: 30
      prefix: "samplerate="
      separate: false

  - id: reads
    type: int?
    inputBinding:
      position: 32
      prefix: "reads="
      separate: false

  - id: mlf
    type: int?
    inputBinding:
      position: 34
      prefix: "mlf="
      separate: false

  - id: minlen
    type: int?
    inputBinding:
      position: 36
      prefix: "minlen="
      separate: false

  - id: restrictleft
    type: int?
    inputBinding:
      position: 38
      prefix: "restrictleft="
      separate: false

  - id: maq
    type: int?
    inputBinding:
      position: 40
      prefix: "maq="
      separate: false

  - id: trimq
    type: int?
    inputBinding:
      position: 42
      prefix: "trimq="
      separate: false

  - id: qtrim
    type: bbduk_types.yaml#side_set?
    inputBinding:
      position: 44
      prefix: "qtrim="
      separate: false

  - id: ftm
    type: int?
    inputBinding:
      position: 46
      prefix: "ftm="
      separate: false

  - id: deterministic
    type: boolean
    default: false
    doc: Produce deterministic output every time
    inputBinding:
      valueFrom: $(bbtoolsSortedOutputsCond(self,inputs.SampleID+inputs.prefix,runtime.cores))
      shellQuote: false
      position: 120

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
    type: File
    outputBinding:
      glob: '*.bbduk.stats'

arguments:

  - valueFrom: $("stats="+inputs.SampleID+inputs.prefix+".bbduk.stats")
  - valueFrom: $("-Xmx"+runtime.ram+"m")

baseCommand: ["bbduk.sh", "statscolumns=5", "ordered"]
