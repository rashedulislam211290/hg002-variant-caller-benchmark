#!/bin/bash
# Step 2: download HG002 chr20 reads and the GIAB v4.2.1 truth set.
# Run from the repository root on a node with internet access (usually the login node):
#   bash 02_download_data.sh
set -eo pipefail
mkdir -p input benchmark

# HG002: Illumina NovaSeq, PCR-free, ~35x, aligned to GRCh38 no-alt, chr20 only
BAM_URL=https://storage.googleapis.com/deepvariant/case-study-testdata
BAM=HG002.novaseq.pcr-free.35x.dedup.grch38_no_alt.chr20.bam
wget -c -P input "${BAM_URL}/${BAM}"
wget -c -P input "${BAM_URL}/${BAM}.bai"

# GIAB HG002 v4.2.1 small-variant benchmark, GRCh38
TRUTH_URL=https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/AshkenazimTrio/HG002_NA24385_son/NISTv4.2.1/GRCh38
wget -c -P benchmark "${TRUTH_URL}/HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz"
wget -c -P benchmark "${TRUTH_URL}/HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz.tbi"
wget -c -P benchmark "${TRUTH_URL}/HG002_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.bed"

samtools quickcheck "input/${BAM}" && echo "BAM OK"
ls -lh input benchmark
