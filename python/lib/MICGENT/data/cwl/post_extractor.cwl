#!/usr/bin/env cwl-runner

cwlVersion: v1.0

class: CommandLineTool
baseCommand: ['python', '-m', 'MICGENT.post_ariba', 'post-extractor','--sample-tars','sample_tars.txt']
requirements:
  - class: InlineJavascriptRequirement
    expressionLib:
      ## this has to be in a function because having it inside InitialWorkDirRequirement
      ## listing:entry causes tabs and newlines quoted due to YAML multiline rules
      - var makeReportList = function(inputs) {
          var s = '';
          for(var i=0; i < inputs.sample_tars.length; ++i) {
              s += inputs.sample_tars[i].path + '\n';
          }
          return s;
        };

  - class: InitialWorkDirRequirement
    listing:
      - entryname: sample_tars.txt
        entry: $(makeReportList(inputs))
doc: |
  Generate post-extraction Web report

inputs:
  manifest:
    type: File
    inputBinding:
        prefix: --manifest
  sequences:
    type: File
    inputBinding:
        prefix: --contigs
  manifest_all:
    type: File
    inputBinding:
        prefix: --manifest-all
  sequences_all:
    type: File
    inputBinding:
        prefix: --contigs-all
  manifest_sum:
    type: File
    inputBinding:
        prefix: --manifest-sum
  manifest_dict:
    type: File
    inputBinding:
        prefix: --manifest-dict
  qc_report:
    type: File
    inputBinding:
        prefix: --qc-report
  qc_data:
    type: File?
    inputBinding:
        prefix: --qc-data
  ref_common:
    type: File
    inputBinding:
        prefix: --ref-common
  micgentjs_tgz:
    type: File
    inputBinding:
        prefix: --micgentjs-tgz
  sor_sfx_target:
    type: string?
    inputBinding:
        prefix: --sor-sfx-target
    default: "X"
  sor_sfx_version:
    type: int?
    inputBinding:
        prefix: --sor-sfx-version
    default: 1
  wf_inputs:
    type: File
    inputBinding:
        prefix: --wf-inputs
  prepareref_tgz:
    type: File
    inputBinding:
        prefix: --prepareref-tgz

  sample_tars: File[]

outputs:
  out_tar:
    type: File
    outputBinding:
      glob: web.tar
  out_sor:
    type: File
    outputBinding:
      glob: sor_pack.tgz
