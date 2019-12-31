class: Workflow
cwlVersion: v1.0
inputs:
  - id: config
    type: File
  - id: sample
    type: replicated_reads.yaml#sample_replicated_reads
  - id: home
    type: string
  - id: threads
    type: int?
  - id: bds_home
    type: string
outputs:
  - id: out_json
    outputSource:
      - atac_encode/out_json
    type: File
  - id: out_archive
    outputSource:
      - atac_encode/out_archive
    type: File
  - id: SampleID
    outputSource:
      - atac_encode/out_SampleID
    type: string
steps:
  - id: cat_files1_1
    in:
      - id: files
        source:
          - sample
        valueFrom: $(self.files1_1)
      - id: out_nameroot
        source:
          - sample
        valueFrom: $(self.SampleID + "R1_1.fastq")
    out:
      - id: output
    run: >-
      cat_files.cwl
  - id: cat_files1_2
    in:
      - id: files
        source:
          - sample
        valueFrom: $(self.files1_2)
      - id: out_nameroot
        source:
          - sample
        valueFrom: $(self.SampleID + "R1_2.fastq")
    out:
      - id: output
    run: >-
      cat_files.cwl
  - id: cat_files2_1
    in:
      - id: files
        source:
          - sample
        valueFrom: $(self.files2_1)
      - id: out_nameroot
        source:
          - sample
        valueFrom: $(self.SampleID + "R2_1.fastq")
    out:
      - id: output
    run: >-
      cat_files.cwl
  - id: cat_files2_2
    in:
      - id: files
        source:
          - sample
        valueFrom: $(self.files2_2)
      - id: out_nameroot
        source:
          - sample
        valueFrom: $(self.SampleID + "R2_2.fastq")
    out:
      - id: output
    run: >-
      cat_files.cwl
  - id: atac_encode
    in:
      - id: SampleID
        source:
          - sample
        valueFrom: $(self.SampleID)
      - id: home
        source:
          - home
      - id: config
        source:
          - config
      - id: threads
        source:
          - threads
      - id: fastq1_1
        source:
          - cat_files1_1/output
      - id: fastq1_2
        source:
          - cat_files1_2/output
      - id: fastq2_1
        source:
          - cat_files2_1/output
      - id: fastq2_2
        source:
          - cat_files2_2/output
      - id: auto_detect_adapter
        default: true
      - id: enable_idr
        default: true
      - id: bds_home
        source:
          - bds_home
    out:
      - id: out_json
      - id: out_archive
      - id: out_SampleID
    run: atac_encode.cwl

requirements:
  - class: StepInputExpressionRequirement
  - class: MultipleInputFeatureRequirement
  - class: SubworkflowFeatureRequirement
  - class: InlineJavascriptRequirement
  - class: SchemaDefRequirement
    types:
      - $import: replicated_reads.yaml
