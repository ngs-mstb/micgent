#!/usr/bin/env cwl-runner

cwlVersion: v1.0

class: CommandLineTool
baseCommand: ['python', '-m', 'MICGENT.gene_extractor', 'combine-samples', 'sample_results.tsv']
requirements:
  - class: InlineJavascriptRequirement
    expressionLib:
      ## this has to be in a function because having it inside InitialWorkDirRequirement
      ## listing:entry causes tabs and newlines quoted due to YAML multiline rules
      - var makeReportList = function(inputs) {
          var s = 'report\tsequence\tstatus\n';
             for(var i=0; i < inputs.reports.length; ++i) {
                s += inputs.reports[i].path + '\t' + inputs.sequences[i].path + 
                '\t' + inputs.statuses[i].path + '\n';
          }
          return s;        
        };
  - class: InitialWorkDirRequirement
    listing:
      - entryname: sample_results.tsv
        entry: $(makeReportList(inputs))
doc: |
  Combine extraction results from multiple samples
arguments:
  - valueFrom: $("manifest_out.tsv")
    position: 101
  - valueFrom: $("seq_out.fasta")
    position: 102


inputs:
  manifest:
    type: File
    doc: Manifest file
    inputBinding:
      position: 1
  reports: File[]
  sequences: File[]
  statuses: File[]

outputs:
  manifest_out:
    type: File
    outputBinding:
      glob: manifest_out.tsv
  sequences_out:
    type: File
    outputBinding:
      glob: seq_out.fasta
