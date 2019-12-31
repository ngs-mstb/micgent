#!/usr/bin/env cwl-runner

class: CommandLineTool
cwlVersion: v1.0

doc: Call our Python all-in-one read trimmin/cleaning/qc-subsampling step

hints:
- class: SoftwareRequirement
  packages:
  - package: 'ngs-mstb'
    version:
    - '1.0'

requirements:
- class: InlineJavascriptRequirement
  expressionLib:
    $import: micgent_js_lib.yaml
- class: SchemaDefRequirement
  types:
    $import: bbduk_types.yaml

inputs:
  - id: SampleID
    type: string

  - id: prefix_out
    type: string?
    default: "_cleaned"

  - id: threads
    type: int?
    default: 1
    inputBinding:
      position: 16
      prefix: "--threads"
      valueFrom: $(runtime.cores)

  - id: ram
    type: int?
    inputBinding:
      position: 16
      prefix: "--ram"
      valueFrom: $(runtime.ram)

  - id: inp_seq1
    type: File
    inputBinding:
      position: 18
      prefix: "--inp-reads"

  - id: inp_seq2
    type: File?
    inputBinding:
      position: 20
      prefix: "--inp-reads2"

  - id: primer_literals
    type: string[]?
    inputBinding:
      position: 22
      prefix: "--primer-literals"
      itemSeparator: ","

  - id: adapter_file
    type: File?
    inputBinding:
      position: 24
      prefix: "--adapter-file"

  - id: clumpify
    type: boolean?
    inputBinding:
      position: 24
      prefix: "--clumpify"

  - id: filter_spikes
    type: boolean?
    inputBinding:
      position: 24
      prefix: "--filter-spikes"

  - id: spikes_file
    type: File?
    inputBinding:
      position: 24
      prefix: "--spikes-file"

  - id: qc_samplerate
    type: float?
    default: 0.1
    inputBinding:
      position: 30
      prefix: "--out-qc-samplerate"

  - id: minlen
    type: int?
    inputBinding:
      position: 36
      prefix: "--minlen"

  - id: maq
    type: int?
    inputBinding:
      position: 40
      prefix: "--maq"

  - id: trimq
    type: int?
    inputBinding:
      position: 42
      prefix: "--trimq"

  - id: qtrim
    type: bbduk_types.yaml#side_set?
    inputBinding:
      position: 44
      prefix: "--qtrim"

  - id: deterministic
    type: boolean
    default: false
    doc: Produce deterministic output every time
    inputBinding:
      prefix: "--deterministic"
      position: 45

  - id: out_seq1
    type: string?
    inputBinding:
      position: 46
      prefix: "--out-reads"
      valueFrom: $( inputs.SampleID+inputs.prefix_out+"_1.fastq.gz" )

  - id: out_seq2
    type: string?
    inputBinding:
      position: 47
      prefix: "--out-reads2"
      valueFrom: >
        ${ if(inputs.inp_seq2) 
          { return inputs.SampleID+inputs.prefix_out+"_2.fastq.gz"; }
          else
          { return null; }
        }

  - id: out_qc_before_seq1
    type: string?
    inputBinding:
      position: 46
      prefix: "--out-qc-before-reads"
      valueFrom: $( inputs.SampleID+inputs.prefix_out+"_1.qc.before.fastq.gz" )

  - id: out_qc_after_seq1
    type: string?
    inputBinding:
      position: 46
      prefix: "--out-qc-after-reads"
      valueFrom: $( inputs.SampleID+inputs.prefix_out+"_1.qc.after.fastq.gz" )

  - id: out_qc_before_seq2
    type: string?
    inputBinding:
      position: 47
      prefix: "--out-qc-before-reads2"
      valueFrom: >
        ${ if(inputs.inp_seq2) 
          { return inputs.SampleID+inputs.prefix_out+"_2.qc.before.fastq.gz"; }
          else
          { return null; }
        }

  - id: out_qc_after_seq2
    type: string?
    inputBinding:
      position: 47
      prefix: "--out-qc-after-reads2"
      valueFrom: >
        ${ if(inputs.inp_seq2) 
          { return inputs.SampleID+inputs.prefix_out+"_2.qc.after.fastq.gz"; }
          else
          { return null; }
        }

  - id: out_stats_base
    type: string?
    inputBinding:
      position: 46
      prefix: "--out-stats"
      valueFrom: $( inputs.SampleID+inputs.prefix_out+".stats" )


outputs:

  - id: out_seqs
    type: File[]
    outputBinding:
      glob: '$(inputs.SampleID)$(inputs.prefix_out)_[0-9].fastq.gz'
      outputEval: $(sortFileObjectsByName(self))

  - id: out_qc_before_seqs
    type: File[]
    outputBinding:
      glob: '$(inputs.SampleID)$(inputs.prefix_out)_[0-9].qc.before.fastq.gz'
      outputEval: $(sortFileObjectsByName(self))

  - id: out_qc_after_seqs
    type: File[]
    outputBinding:
      glob: '$(inputs.SampleID)$(inputs.prefix_out)_[0-9].qc.after.fastq.gz'
      outputEval: $(sortFileObjectsByName(self))

  - id: out_SampleID
    type: string
    outputBinding:
      outputEval: $(inputs.SampleID)

  - id: out_stats
    type: File
    outputBinding:
      glob: $( inputs.SampleID+inputs.prefix_out+".stats" )

arguments:

baseCommand: ["python","-m","MICGENT.pysteps","clean-reads","--stdout","stdout","--stderr","stderr"]
