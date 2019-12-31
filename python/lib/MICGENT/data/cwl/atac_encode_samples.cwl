class: Workflow
cwlVersion: v1.0
inputs:
  - id: samples
    type: replicated_reads.yaml#sample_replicated_reads[]
  - id: atac_encode__home
    type: string
  - id: atac_encode__threads
    type: int?
  - id: bds_home
    type: string
  - id: atac_encode__config
    type: File

outputs:
  - id: out_jsons
    outputSource:
      - atac_encode_sample/out_json
    type: File[]
  - id: out_archives
    outputSource:
      - atac_encode_sample/out_archive
    type: File[]
  - id: SampleIDs
    outputSource:
      - atac_encode_sample/SampleID
    type: string[]
steps:
  - id: atac_encode_sample
    in:
      - id: sample
        source:
          - samples
      - id: home
        source:
          - atac_encode__home
      - id: bds_home
        source:
          - bds_home
      - id: config
        source:
          - atac_encode__config
      - id: threads
        source:
          - atac_encode__threads
    out:
      - id: out_json
      - id: out_archive
      - id: SampleID
    scatter:
      - sample
    run: >-
      atac_encode_sample.cwl
requirements:
  - class: StepInputExpressionRequirement
  - class: MultipleInputFeatureRequirement
  - class: ScatterFeatureRequirement
  - class: SubworkflowFeatureRequirement
  - class: InlineJavascriptRequirement
  - class: SchemaDefRequirement
    types:
      - $import: replicated_reads.yaml
