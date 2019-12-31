class: Workflow
cwlVersion: v1.0
hints:
  SoftwareRequirement:
    packages:
    - package: 'ngs-mstb'
      version:
      - '1.0'
inputs:
  - id: prepareref_tgz
    type: File
  - id: manifest
    type: File
  - id: samples
    type: 'sample_reads.yaml#sample_reads[]'
  - id: ariba_threads
    type: int?
    ## Note that we need to compute ariba_threads in workflow step input,
    ## because here in workflow inputs valueFrom is slidently ignored by
    ## both cwltool and cwltoil.
  - id: trim_threads
    type: int?
    default: 4
  - id: qc_threads
    type: int?
  - id: qc_samplerate
    type: float?
    default: 0.1
  - id: clumpify
    type: boolean?
  - id: filter_spikes
    type: boolean?
  - id: primer_literals
    type: string[]?
  - id: adapter_file
    type: File?
  - id: read_minlen
    type: int?
    default: 50
  - id: spikes_file
    type: File?
  - id: assembly_cov
    type: int?
  - id: assembly_cov_min
    type: int?
  - id: gene_nt_extend
    type: int?
  - id: nucmer_min_id
    type: int?
  - id: nucmer_min_len
    type: int?
  - id: nucmer_breaklen
    type: int?
  - id: assembler
    type: string?
    default: "fermilite"
  - id: spades_mode
    type: string?
  - id: spades_options
    type: string?
  - id: plugin_asm_options
    type: string?
  - id: ref_common
    type: File
  - id: micgentjs_tgz
    type: File
  - id: filter_asm_args
    type: string?
  - id: sor_sfx_target
    type: string?
    default: "X"
  - id: sor_sfx_version
    type: int?
    default: 1
  - id: wf_inputs
    type: File
  - id: debug
    type: boolean
    default: false
    ## this will cause various efficiency losses, so the default is false
  - id: deterministic
    type:  boolean?
    default: false
  - id: cut_to_ref
    type: boolean?
    default: false
  - id: one_seq_per_contig
    type: boolean?
    default: false
  - id: pad_assembled
    type: int?
    default: 200
  - id: pad_gene
    type: int?
    default: 200
  - id: sig_inp_key
    type: string?
    default: "123"

outputs:
  - id: manifest_out
    outputSource:
      - filter_assemblies/manifest_out
    type: File
  - id: sequences_out
    outputSource:
      - filter_assemblies/sequences_out
    type: File
  - id: clean_stats
    outputSource:
      - all_samples/clean_stats
    type: File[]
  - id: qc_report
    outputSource:
      - multiqc/out_report
    type: File
  - id: qc_data
    outputSource:
      - multiqc/out_data
    type: File
  - id: cleaned_seqs
    outputSource:
      - flatten_cleaned_reads/flat
    type: File[]
  - id: web_tar_out
    outputSource:
      - post_extractor/out_tar
    type: File
  - id: sor_pack_out
    outputSource:
      - post_extractor/out_sor
    type: File
  - id: ariba_debug_tar
    outputSource:
      - all_samples/ariba_debug_tar
    type: File[]?


steps:
  - id: all_samples
    in:
      - id: prepareref_tgz
        source:
          - prepareref_tgz
      - id: sample
        source:
          - samples
      - id: ariba_threads
        source:
          - ariba_threads
      - id: trim_threads
        source:
          - trim_threads
      - id: qc_threads
        source:
          - qc_threads
      - id: qc_samplerate
        source:
          - qc_samplerate
      - id: clumpify
        source:
          - clumpify
      - id: filter_spikes
        source:
          - filter_spikes
      - id: primer_literals
        source:
          - primer_literals
      - id: adapter_file
        source:
          - adapter_file
      - id: spikes_file
        source:
          - spikes_file
      - id: read_minlen
        source:
          - read_minlen
      - id: assembly_cov
        source:
          - assembly_cov
      - id: assembly_cov_min
        source:
          - assembly_cov_min
      - id: gene_nt_extend
        source:
          - gene_nt_extend
      - id: nucmer_min_id
        source:
          - nucmer_min_id
      - id: nucmer_min_len
        source:
          - nucmer_min_len
      - id: nucmer_breaklen
        source:
          - nucmer_breaklen
      - id: assembler
        source:
          - assembler
      - id: spades_mode
        source:
          - spades_mode
      - id: spades_options
        source:
          - spades_options
      - id: plugin_asm_options
        source:
          - plugin_asm_options
      - id: debug
        source:
          - debug
      - id: deterministic
        source:
          - deterministic
      - id: cut_to_ref
        source:
          - cut_to_ref
      - id: one_seq_per_contig
        source:
          - one_seq_per_contig
      - id: pad_assembled
        source:
          - pad_assembled
      - id: pad_gene
        source:
          - pad_gene
      - id: sig_inp_key
        source:
          - sig_inp_key

    out: [report, sequences, status, clean_stats, qc_reports_before, qc_reports_after, cleaned_seqs, asmqc_tar, ariba_debug_tar]
    run:
      class: Workflow
      cwlVersion: v1.0
      inputs:
        - id: prepareref_tgz
          type: File
        - id: sample
          type: 'sample_reads.yaml#sample_reads'
        - id: ariba_threads
          type: int?
        - id: trim_threads
          type: int?
        - id: qc_threads
          type: int?
        - id: qc_samplerate
          type: float?
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
        - id: read_minlen
          type: int?
        - id: assembly_cov
          type: int?
        - id: assembly_cov_min
          type: int?
        - id: gene_nt_extend
          type: int?
        - id: nucmer_min_id
          type: int?
        - id: nucmer_min_len
          type: int?
        - id: nucmer_breaklen
          type: int?
        - id: assembler
          type: string?
        - id: spades_mode
          type: string?
        - id: spades_options
          type: string?
        - id: plugin_asm_options
          type: string?
        - id: debug
          type: boolean
        - id: deterministic
          type:  boolean?
          default: false
        - id: cut_to_ref
          type: boolean?
          default: false
        - id: one_seq_per_contig
          type: boolean?
          default: false
        - id: pad_assembled
          type: int?
          default: 200
        - id: pad_gene
          type: int?
          default: 200
        - id: sig_inp_key
          type: string

      outputs:
        - id: report
          outputSource:
            - extract_contigs/report_out
          type: File
        - id: sequences
          outputSource:
            - extract_contigs/sequences
          type: File
        - id: status
          outputSource:
            - extract_contigs/status_out
          type: File
        - id: clean_stats
          outputSource:
            - clean_reads/out_stats
          type: File
        - id: qc_reports_before
          outputSource:
            - clean_reads/out_qc_reports_before
          type: File[]
        - id: qc_reports_after
          outputSource:
            - clean_reads/out_qc_reports_after
          type: File[]
        - id: cleaned_seqs
          outputSource:
            - clean_reads/out_seqs
          type: File[]
        - id: asmqc_tar
          outputSource:
            - ariba_run/asmqc_tar
          type: File
        - id: ariba_debug_tar
          outputSource:
            - ariba_run/ariba_debug_tar
          type: File?

      steps:

        - id: sig_inp1
          in:
            inp_file:
              source: sample
              valueFrom: $(self.file1)
            out_copy_prefix:
              valueFrom: "imp_"
            key: sig_inp_key
          out:
            - sig
            - out_copy
          run: sig.cwl
          requirements:
            - class: ResourceRequirement
              coresMax: 1

        - id: sig_inp2
          in:
            inp_file:
              source: sample
              valueFrom: $(self.file2)
            out_copy_prefix:
              valueFrom: "imp_"
            key: sig_inp_key
          out:
            - sig
            - out_copy
          run: sig.cwl
          requirements:
            - class: ResourceRequirement
              coresMax: 1

        - id: clean_reads
          in:
            SampleID:
              source: sample
              valueFrom: $(self.SampleID)
            inp_seq1: sig_inp1/out_copy
            inp_seq2: sig_inp2/out_copy
            qc_samplerate: qc_samplerate
            clumpify: clumpify
            filter_spikes: filter_spikes
            primer_literals: primer_literals
            adapter_file: adapter_file
            spikes_file: spikes_file
            minlen: read_minlen
            threads: trim_threads
            qc_threads: qc_threads
            deterministic: deterministic
            prefix_out:
              valueFrom: "_cleaned"
          out:
            - out_seqs
            - out_stats
            - out_qc_reports_before
            - out_qc_reports_after
          requirements:
            - class: ResourceRequirement
              coresMin: $(inputs.threads)
              ramMin: 4096
          run: clean_reads_qc.cwl

        - id: ariba_run
          in:
            SampleID:
              source: sample
              valueFrom: $(self.SampleID)
            prepareref_tgz: prepareref_tgz
            reads_1:
              source: clean_reads/out_seqs
              valueFrom: $(self[0])
            reads_2:
              source: clean_reads/out_seqs
              valueFrom: $(self[1])
            assembly_cov: assembly_cov
            assembly_cov_min: assembly_cov_min
            assembler: assembler
            spades_mode: spades_mode
            spades_options: spades_options
            plugin_asm_options: plugin_asm_options
            threads:
              source: ariba_threads
              valueFrom: $(getAribaCores(self,runtime,inputs))
            debug: debug
            serial: deterministic
            gene_nt_extend: gene_nt_extend
            nucmer_min_id: nucmer_min_id
            nucmer_min_len: nucmer_min_len
            nucmer_breaklen: nucmer_breaklen

          out:
            - assembled_genes
            - assembled_seqs
            - assemblies
            - log_clusters
            - report
            - version_info
            - asmqc_tar
            - ariba_debug_tar
            - seq_map
            - basecov_asm
            - basecov_ref
            - stdout
            - stderr
          run: ariba_run.cwl
          requirements:
            - class: ResourceRequirement
              coresMin: $(inputs.threads)
              ramMin: 8192
        - id: extract_contigs
          in:
            cut_to_ref: cut_to_ref
            one_seq_per_contig: one_seq_per_contig
            pad_assembled: pad_assembled
            pad_gene: pad_gene
            SampleID:
              source: sample
              valueFrom: $(self.SampleID)
            report: ariba_run/report
            assembled_genes: ariba_run/assembled_genes
            assembled_seqs: ariba_run/assembled_seqs
            assemblies: ariba_run/assemblies
            seq_map: ariba_run/seq_map
            ariba_stdout: ariba_run/stdout
            ariba_stderr: ariba_run/stderr
            basecov_asm: ariba_run/basecov_asm
            basecov_ref: ariba_run/basecov_ref
            sig_inp: [ sig_inp1/sig, sig_inp2/sig ]
          out:
            - report_out
            - sequences
            - status_out
          run: extract_contigs.cwl
          requirements:
            - class: ResourceRequirement
              coresMax: 1
      requirements:
        - class: InlineJavascriptRequirement
          expressionLib:
            - $import: micgent_js_lib.yaml        
        - class: StepInputExpressionRequirement
        - class: MultipleInputFeatureRequirement

    scatter:
      - sample

  - id: combine_samples
    in:
      manifest: manifest
      reports: all_samples/report
      sequences: all_samples/sequences
      statuses: all_samples/status
    out:
      - manifest_out
      - sequences_out
    run: combine_samples.cwl
    requirements:
      - class: ResourceRequirement
        coresMax: 1

  - id: filter_assemblies
    in:
      manifest: combine_samples/manifest_out
      sequences: combine_samples/sequences_out
      args: filter_asm_args
    out:
      - manifest_out
      - sequences_out
      - manifest_out_all
      - manifest_out_sum
      - manifest_out_dict
    run: filter_assemblies.cwl
    requirements:
      - class: ResourceRequirement
        coresMax: 1

  - id: flatten_qc
    in:
      inp:
        source: [all_samples/qc_reports_before, all_samples/qc_reports_after]
        linkMerge: merge_flattened
    out: [ flat ]
    run: et_flatten_array.cwl

  - id: flatten_cleaned_reads
    in:
      inp: all_samples/cleaned_seqs
    out: [ flat ]
    run: et_flatten_array.cwl

  - id: multiqc
    in:
      inp_files: flatten_qc/flat
      config_str:
        valueFrom: >
          module_order:
            - fastqc:
                name: 'FastQC before read cleaning'
                path_filters:
                  - '*before_cleaning_*_fastqc.zip'
            - fastqc:
                name: 'FastQC after read cleaning'
                path_filters:
                  - '*after_cleaning_*_fastqc.zip'
      filename:
        valueFrom: "multiqc"
    out: [out_report, out_data]
    run: multiqc.cwl
    requirements:
      - class: ResourceRequirement
        coresMax: 1

  - id: post_extractor
    in:
      manifest: filter_assemblies/manifest_out
      sequences: filter_assemblies/sequences_out
      manifest_all: filter_assemblies/manifest_out_all
      sequences_all: combine_samples/sequences_out
      manifest_sum: filter_assemblies/manifest_out_sum
      manifest_dict: filter_assemblies/manifest_out_dict
      qc_report: multiqc/out_report
      qc_data: multiqc/out_data
      ref_common: ref_common
      micgentjs_tgz: micgentjs_tgz
      sample_tars: all_samples/asmqc_tar
      prepareref_tgz: prepareref_tgz
      wf_inputs: wf_inputs
      sor_sfx_target: sor_sfx_target
      sor_sfx_version: sor_sfx_version

    out: [out_tar,out_sor]
    run: post_extractor.cwl
    requirements:
      - class: ResourceRequirement
        coresMax: 1

requirements:
  - class: ScatterFeatureRequirement
  - class: InlineJavascriptRequirement
    expressionLib:
      - $import: micgent_js_lib.yaml

  - class: StepInputExpressionRequirement
  - class: SubworkflowFeatureRequirement
  - class: MultipleInputFeatureRequirement
  - class: SchemaDefRequirement
    types:
      - $import: sample_reads.yaml
