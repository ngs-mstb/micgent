#!/usr/bin/env cwl-runner

class: CommandLineTool
cwlVersion: v1.0

doc: Use Clumpify from BBMap toolset to remove (optical) duplicates

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
  - id: SampleID
    type: string

  - id: prefix?
    type: string
    default: "_cl"

  - id: dupedist
    type: int?
    inputBinding:
      position: 8
      prefix: "dupedist="
      separate: false

  - id: dedupe
    type:  boolean?
    inputBinding:
      position: 10
      prefix: "dedupe"
    default: true

  - id: optical
    type:  boolean?
    inputBinding:
      position: 12
      prefix: "optical"
    default: true

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
      ## as of version 37.99, the output is non-deterministic in multithreaded mode, and
      ## deterministic with a single thread, but sorting fixes it with multithreads
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
    type: stdout

arguments:

  - valueFrom: $("-Xmx"+runtime.ram+"m")

baseCommand: ["clumpify.sh"]
stdout: $(inputs.SampleID+inputs.prefix+".clumpify.stats")