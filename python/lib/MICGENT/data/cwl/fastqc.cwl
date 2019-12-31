#!/usr/bin/env cwl-runner

cwlVersion: v1.0

doc: Generate report with FASTQC

requirements:
- class: ShellCommandRequirement
- class: InlineJavascriptRequirement
  expressionLib:
    $import: micgent_js_lib.yaml
- class: EnvVarRequirement
  envDef:
    _JAVA_OPTIONS: "-Xms1024m -Xmx$(runtime.ram)m -XX:ParallelGCThreads=1"

class: CommandLineTool

inputs:
  - id: adapters
    type: File?
    inputBinding:
      prefix: --adapters

  ## What we need for the *_fastqc.html files here and in the output glob below:
  ## If given input with zero reads, FASTQ Java program will traceback, but Perl wrapper exit with code 0,
  ## and a corrupted zip file will be created. But the HTML file will not be created. So, we use the presence of
  ## HTML file as a flag that the corresponding .zip file should be returned or not.
  - id: prefix
    type: string?
    inputBinding:
      valueFrom: ${ if(self) { return '; prefix="'+self+'"; (for f in *_fastqc.zip *_fastqc.html; do [ -f "$f" ] && mv "$f" "${prefix}$f"; done); true '; } else { return '; true '; } }
      shellQuote: false
      position: 1000

  - id: casava
    type: boolean
    default: false
    inputBinding:
      prefix: --casava

  - id: contaminants
    type: File?
    inputBinding:
      prefix: --contaminants

  - id: extract
    type: boolean
    default: false
    inputBinding:
      prefix: --extract

  - id: format
    type: string
    default: fastq
    inputBinding:
      prefix: --format

  - id: inp_seqs
    type: File[]
    inputBinding:
      position: 99

  - id: kmers
    type: File?
    inputBinding:
      prefix: --kmers

  - id: limits
    type: File?
    inputBinding:
      prefix: --limits

  - id: nano
    type: boolean
    default: false
    inputBinding:
      prefix: --nano

  - id: noextract
    type: boolean
    default: true
    inputBinding:
      prefix: --noextract

  - id: nofilter
    type: boolean
    default: false
    inputBinding:
      prefix: --nofilter

  - id: nogroup
    type: boolean
    default: false
    inputBinding:
      prefix: --nogroup

  - id: quiet
    type: boolean
    default: false
    inputBinding:
      prefix: --quiet

  - id: threads
    type: int?
    inputBinding:
      prefix: --threads
      valueFrom: $(getCores(self,runtime))

outputs:
  - id: out_report
    type: File[]
    outputBinding:
      glob: "*_fastqc.html"
      outputEval: $(self.map(x => ({"class":"File","path":x.nameroot+'.zip'})))

arguments:
  - prefix: --dir
    valueFrom: $(runtime.tmpdir)
  - prefix: --outdir
    valueFrom: $(runtime.outdir)

baseCommand: ["fastqc"]
