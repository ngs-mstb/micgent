#!/usr/bin/env cwl-runner

class: CommandLineTool
cwlVersion: v1.0

doc: Use BBmap from BBTools to map reads

hints:
- class: ResourceRequirement
  ramMin: 16384
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
    default: "_mp"

  - id: nodisk
    type:  boolean?
    inputBinding:
      position: 10
      prefix: "nodisk"

  - id: ref
    type: File?
    inputBinding:
      position: 18
      prefix: "ref="
      separate: false

  - id: ref_ind
    #This would better be defined with type: Directory, but Toil does not support Directory type,
    #so this has to be a string containing a path to a directory that is local to
    #the execution node, and pre-built outside of CWL. These directories are not small
    #(4GB per human reference), so making them into tarballs that are pushed into the
    #execution as File and unpacked would be too expensive. Tarball might ge cached, but not the
    #unpacked version.
    #Note: secondaryFile workaround will not work for this because bbmap's ref has deep structure.
    #Toil supports InitialWorkDirRequirement and ExpressionTool with some limitations.
    # Arranging files with symlinks as part of tool
    #execution will probably work. It should be possible to write a universal Python script that
    #creates a flattened symlinks after the main command, and unflattens secondaryFiles before the
    #main command. Maybe even better to push this into Toil so that it finally supports the Directory
    #requirement. BCBio does something along those lines. Actually, do this:
    #Output record that contains String[],File[], and has outputBinding with a glob with nested wildcards (dir/*/*/*)
    #to build the output array, and the String element has outputEval to strip the leading directory component.
    #On input, InitialWorkDirRequirement (or CreateFileRequirement) can simply generate a shell script
    #with as many lines as there are files that calls mkdir -p and ln -s. That script is then executed
    #as the first part of the tool command.
    type: string?
    inputBinding:
      position: 18
      prefix: "path="
      separate: false

  - id: maxindel
    type: int?
    inputBinding:
      position: 8
      prefix: "maxindel="
      separate: false

  - id: minid
    type: float?
    inputBinding:
      position: 8
      prefix: "minid="
      separate: false

  - id: bwr
    type: float?
    inputBinding:
      position: 8
      prefix: "bwr="
      separate: false

  - id: bw
    type: int?
    inputBinding:
      position: 8
      prefix: "bw="
      separate: false

  - id: minhits
    type: int?
    inputBinding:
      position: 8
      prefix: "minhits="
      separate: false

  - id: quickmatch
    type:  boolean?
    inputBinding:
      position: 10
      prefix: "quickmatch"

  - id: untrim
    type:  boolean?
    inputBinding:
      position: 10
      prefix: "untrim"

  - id: fast
    type:  boolean?
    inputBinding:
      position: 10
      prefix: "fast"

  - id: trimq
    type: int?
    inputBinding:
      position: 12
      prefix: "trimq="
      separate: false

  - id: qtrim
    type: bbduk_types.yaml#side_set?
    inputBinding:
      position: 12
      prefix: "qtrim="
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

  - id: outu
    type: string?
    inputBinding:
      position: 20
      prefix: "outu="
      separate: false
      valueFrom: $(inputs.SampleID+inputs.prefix+"_"+self)

  - id: outu2
    type: string?
    inputBinding:
      position: 20
      prefix: "outu2="
      separate: false
      valueFrom: $(inputs.SampleID+inputs.prefix+"_"+self)

  - id: outm
    type: string?
    inputBinding:
      position: 20
      prefix: "outm="
      separate: false
      valueFrom: $(inputs.SampleID+inputs.prefix+"_"+self)

  - id: outm2
    type: string?
    inputBinding:
      position: 20
      prefix: "outm2="
      separate: false
      valueFrom: $(inputs.SampleID+inputs.prefix+"_"+self)


outputs:

  - id: outs_u
    type: File[]?
    outputBinding:
      glob: ['*$(inputs.outu)','*$(inputs.outu2)']
      outputEval: $(sortFileObjectsByName(self))

  - id: outs_m
    type: File[]?
    outputBinding:
      glob: ['*$(inputs.outm)','*$(inputs.outm2)']
      outputEval: $(sortFileObjectsByName(self))

  - id: out_SampleID
    type: string
    outputBinding:
      outputEval: $(inputs.SampleID)

  - id: out_stats
    type: stdout

arguments:

  - valueFrom: $("-Xmx"+runtime.ram+"m")

baseCommand: ["bbmap.sh"]
stdout: $(inputs.SampleID+inputs.prefix+".bbmap.stats")