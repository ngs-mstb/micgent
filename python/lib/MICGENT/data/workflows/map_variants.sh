#!/bin/bash
set -ex
[ -n "$1" ] || exit 1
ref=$1
shift
[ -n "$1" ] || exit 1
bam=$1
shift
[ -n "$1" ] || exit 1
bcf=$1
shift
[ -n "$1" ] || exit 1
threads=$1
shift
[ -n "$1" ] || exit 1
read1=$1
shift
[ -n "$1" ] || exit 1
read2=$1

tmp_root=map
ln -sf "$ref" "${tmp_root}.ref.fa"
ref="${tmp_root}.ref.fa"

sickle pe -f $read1 -r $read2 -t sanger -o ${tmp_root}.1.fq -p ${tmp_root}.2.fq -s ${tmp_root}.s.fq -q 20 -l 50
read1=${tmp_root}.1.fq
read2=${tmp_root}.2.fq
reads="$read1 $read2"

bwa index $ref
bwa mem -t "$threads" -M "$ref" $reads \
    | samtools view -@ "$threads" -bSu - \
    | samtools sort -@ "$threads" -T ${tmp_root}.tmp -o ${tmp_root}.bam -
samtools index ${tmp_root}.bam
picard -XX:ParallelGCThreads=1 MarkDuplicates \
                            MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=1000 \
                            METRICS_FILE=${tmp_root}.dedup.metrics \
                            REMOVE_DUPLICATES=true \
                            ASSUME_SORTED=true  \
                            VALIDATION_STRINGENCY=LENIENT \
                            INPUT=${tmp_root}.bam \
                            OUTPUT="$bam"
samtools index "$bam"

#let "threads_bcf=$threads-1"
threads_bcf=$threads
samtools mpileup -uf "$ref" "$bam" \
    | bcftools call --threads "$threads_bcf" -vm - > "$bcf"

