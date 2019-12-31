#!/bin/sh
this_dir=$(cd $(dirname $0); pwd)
deploy_root="$1"
shift
[ -n "$deploy_root" ] || exit 1
pytest --run-slow \
    --ignore test_run \
    --conda-env-ngs-mstb=ngs-mstb \
    --micgent-data=$deploy_root/micgent_db \
    --extra-config=$this_dir/run_tests_gene_extractor.yaml \
    --log-cli-level=DEBUG \
    --log-file=test_std.log \
    --large-test-data=$deploy_root/micgent_large_test_data
#    -k test_run_extractor_determinism_rsv_1
#    -k test_run_extractor_small_rsv
#    --huge-test-data=$deploy_root/micgent_huge_test_data \
#    -k test_run_extractor_determinism_rsv_complete_1
#    -k test_run_extractor_small_sa
#    -k test_filter_assemblies
#    -k test_simplify_rerun_json
#    -k test_clean_reads_small_sa
#    -k test_run_extractor_zero_reads_rsv
#    -k test_extract_contigs_skewed_cov_filter_rsv
#    -k test_extract_contigs
#    -k test_file_sig_copy
#    -k test_file_sig_cmp_msg
#    -k test_run_extractor_chim_duplicate_rsv
#    -k test_run_extractor_no_min_cov_reads_rsv
#    -k test_run_extractor_skewed_cov_filter_rsv
#    -k test_ariba_skewed_cov_filter_rsv
#    -k 'test_run_extractor_determinism_multi_run_huge_sa_2'
#    -k test_run_extractor_zero_reads_rsv
#    -k test_ariba_no_min_cov_reads_rsv
#    -k test_run_extractor_stop_codon_sa
#    -k test_extract_contigs_stop_codon_revcomp
#    -k test_ariba_stop_codon_sa
#    -k test_extract_contigs_stop_codon
#    -k test_ariba_small_sa
#    -k 'test_run_extractor_determinism_multi_run_huge_sa_2'
#    -k 'test_run_extractor_determinism_multi_run_huge_rsv_1'
#    -k test_run_extractor_determinism_multi_run_rsv_1
#    -k test_ariba_zero_reads_rsv
#    -k test_ariba_one_read_rsv
#    -k est_run_extractor_determinism_sa_fermi_2
#    -k test_run_extractor_determinism_sa_2
#    -k test_filter_assemblies
#    -k test_run_toil
#    -k test_run_extractor_determinism_multi_run_sa_2
#    -k test_run_extractor_determinism_sa_2
#    -k test_clean_reads_full_size_sa
#    -k test_clean_reads_small_rsv
#    -k test_clean_reads_rsv
#    -k 'test_run_extractor_determinism_multi_run_huge_'
#    -k test_run_extractor_determinism
#    -k test_gene_extractor
