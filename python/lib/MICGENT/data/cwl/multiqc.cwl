#!/usr/bin/env cwl-runner

cwlVersion: v1.0

doc: Generate combined report with MULTIQC

hints:
- class: ResourceRequirement
  ramMin: 2048

requirements:
- class: InlineJavascriptRequirement
  expressionLib:
    ## this has to be in a function because having it inside InitialWorkDirRequirement
    ## listing:entry causes tabs and newlines quoted due to YAML multiline rules
    - var makeReportList = function(inputs) {
        var out = ['empty.txt'];
        inputs.inp_files.forEach(function(x) {
            out.push(x.path);
        });
        return out.join('\n');
      };
- class: EnvVarRequirement
  envDef:
    LANG: en_US.UTF-8
    LC_ALL: en_US.UTF-8
- class: InitialWorkDirRequirement
  listing:
  - entryname: inp_list.txt
    entry: $(makeReportList(inputs))
  - entryname: config.yaml
    entry: ${ if(inputs.config_str) { return inputs.config_str; } else { return ''; } }
  - entryname: empty.txt
    entry: ''
  - entryname: $(inputs.filename).html
    entry: |
      <html lang="en">
      <body>
      No sample data available to build a report
      </body>
      </html>

class: CommandLineTool

inputs:
  - id: inp_files
    type: File[]

  - id: config
    type: File?

  - id: config_str
    type: string?

  - id: cl_config
    type: string?

  - id: filename
    type: string?
    inputBinding:
      prefix: --filename

  - id: data_dir
    type: boolean
    default: true
    inputBinding:
      prefix: --data-dir

  - id: zip_data_dir
    type: boolean
    default: true
    inputBinding:
      prefix: --zip-data-dir

  - id: interactive
    type: boolean
    default: false
    inputBinding:
      prefix: --interactive

  - id: flat
    type: boolean
    default: false
    inputBinding:
      prefix: --flat

outputs:
  - id: out_report
    type: File
    outputBinding:
      glob: $(inputs.filename).html
  - id: out_data
    type: File?
    outputBinding:
      glob: $(inputs.filename)_data.zip

arguments:
  - prefix: --config
    valueFrom: ${ if(inputs.config_str) { return "config.yaml"; } else if(inputs.config) { return inputs.config.path; } else {return null;} }
    #valueFrom: $({"class":"File","path":"config.yaml"})


baseCommand: ["multiqc","--force","--file-list","inp_list.txt"]
