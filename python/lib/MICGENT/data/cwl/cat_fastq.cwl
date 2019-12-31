class: Workflow
cwlVersion: v1.0
inputs:
  - id: samples
    type: 'multi_reads.yaml#sample_multi_reads[]'
outputs:
  - id: files1
    type: File[]
    outputSource:
      - all_samples/file1
  - id: files2
    type: File[]
    outputSource:
      - all_samples/file2
steps:
  - id: all_samples
    in:
      - id: sample
        source:
          - samples
    out: [file1,file2]
    scatter:
      - sample
    run:
      class: Workflow
      cwlVersion: v1.0
      inputs:
        - id: sample
          type: multi_reads.yaml#sample_multi_reads
      outputs:
        - id: file1
          outputSource:
            - cat_files1/output
          type: File
        - id: file2
          outputSource:
            - cat_files2/output
          type: File
      steps:
        - id: cat_files1
          in:
            - id: files
              source:
                - sample
              valueFrom: $(self.files1)
            - id: out_nameroot
              source:
                - sample
              valueFrom: $(self.SampleID + "_1.fastq")
          out:
            - id: output
          run: >-
            cat_files.cwl
        - id: cat_files2
          in:
            - id: files
              source:
                - sample
              valueFrom: $(self.files2)
            - id: out_nameroot
              source:
                - sample
              valueFrom: $(self.SampleID + "_2.fastq")
          out:
            - id: output
          run: >-
            cat_files.cwl

requirements:
  - class: StepInputExpressionRequirement
  - class: MultipleInputFeatureRequirement
  - class: SubworkflowFeatureRequirement
  - class: InlineJavascriptRequirement
  - class: ScatterFeatureRequirement
  - class: SchemaDefRequirement
    types:
      - $import: multi_reads.yaml
