class: Workflow
cwlVersion: v1.0

inputs:
  - id: SampleID
    type: string

  - id: threads
    type: int?
    default: 4

  - id: qc_threads
    type: int?
    default: 0

  - id: inp_seq1
    type: File

  - id: inp_seq2
    type: File?

  - id: clumpify
    type: boolean?

  - id: filter_spikes
    type: boolean?

  - id: primer_literals
    type: string[]?

  - id: adapter_file
    type: File?

  - id: spikes_file
    type: File?

  - id: qc_samplerate
    type: float?
    default: 0.1

  - id: minlen
    type: int?
    # Before FASTQC v.0.11.7, we were hitting the bug below. It makes sense to drop very short
    # reads anyway after trimming, so we keep the default minlen filter.
    # https://github.com/s-andrews/FastQC/issues/1
    default: 13

  - id: maq
    type: int?

  - id: trimq
    type: int?

  - id: qtrim
    type: bbduk_types.yaml#side_set?

  - id: deterministic
    type:  boolean?
    default: false

  - id: prefix_out
    type: string?
    default: "_cleaned"

outputs:

  - id: out_seqs
    type: File[]
    outputSource:
      clean_reads/out_seqs

  - id: out_SampleID
    type: string
    outputSource:
      clean_reads/out_SampleID

  - id: out_stats
    type: File
    outputSource:
      clean_reads/out_stats

  - id: out_qc_reports_before
    type: File[]
    outputSource:
      qc_before/out_report

  - id: out_qc_reports_after
    type: File[]
    outputSource:
      qc_after/out_report

steps:
  - id: clean_reads
    in:
      - id: inp_seq1
        source:
          - inp_seq1
      - id: inp_seq2
        source:
          - inp_seq2
      - id: SampleID
        source:
          - SampleID
      - id: clumpify
        source: clumpify
      - id: filter_spikes
        source: filter_spikes
      - id: primer_literals
        source: primer_literals
      - id: adapter_file
        source: adapter_file
      - id: qc_samplerate
        source: qc_samplerate
      - id: minlen
        source: minlen
      - id: maq
        source: maq
      - id: trimq
        source: trimq
      - id: qtrim
        source: qtrim
      - id: deterministic
        source:
          - deterministic
      - id: prefix_out
        source:
          - prefix_out
      - id: threads
        source: threads
    out:
      - id: out_seqs
      - id: out_SampleID
      - id: out_stats
      - id: out_qc_before_seqs
      - id: out_qc_after_seqs

    run: >-
      clean_reads.cwl
    requirements:
      - class: ResourceRequirement
        coresMin: $(inputs.threads)
        ramMin: 4096

  - id: qc_before
    in:
      - id: inp_seqs
        source: clean_reads/out_qc_before_seqs
      - id: prefix
        valueFrom: "before_cleaning_"
      ## what this is for: to access both thread arguments in coresMin below, because
      ## 'inputs' in the step relates to step inputs, not workflow inputs. And
      ## valueFrom in workflow inputs is simply ignored.
      - id: threads
        source: [qc_threads,threads]
        valueFrom: $(getFastqcCores(self,runtime))
    out:
      - id: out_report
    run: >-
      fastqc.cwl
    ## TODO: expose those dedup parameters that have to be changed by the user according to the instrument model
    ## https://www.biostars.org/p/225338/
    ## For dedup, as the author advises, it is best to keep subs param at a non-zero value, otherwise will enrich dataset
    ## with reads containing errors.
    requirements:
      - class: ResourceRequirement
        coresMin: $(inputs.threads)
        ramMin: 2048

  - id: qc_after
    in:
      - id: inp_seqs
        source:
          - clean_reads/out_qc_after_seqs
      - id: prefix
        valueFrom: "after_cleaning_"
      - id: threads
        source: [qc_threads,threads]
        valueFrom: $(getFastqcCores(self,runtime))
    out:
      - id: out_report
    run: >-
      fastqc.cwl
    requirements:
      - class: ResourceRequirement
        coresMin: $(inputs.threads)
        ramMin: 2048

requirements:
  - class: StepInputExpressionRequirement
  - class: SubworkflowFeatureRequirement
  - class: MultipleInputFeatureRequirement
  - class: InlineJavascriptRequirement
  - class: SchemaDefRequirement
    types:
      $import: bbduk_types.yaml
  - class: InlineJavascriptRequirement
    expressionLib:
      $import: micgent_js_lib.yaml

hints:
  - class: ResourceRequirement
    coresMin: $(inputs.threads)
  - class: SoftwareRequirement
    packages:
    - package: 'ngs-mstb'
      version:
      - '1.0'
