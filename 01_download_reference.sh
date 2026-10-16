#!/bin/bash
# Step 1: download and index the GRCh38 no-alt analysis set.
# Run from the repository root on a node with internet access (usually the login node):
#   bash scripts/01_download_reference.sh
set -eo pipefail

REF=GCA_000001405.15_GRCh38_no_alt_analysis_set.fasta
URL=https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/references/GRCh38/${REF}.gz

if [ ! -s "${REF}" ]; then
  wget -c "${URL}"
  gunzip "${REF}.gz"
fi
samtools faidx "${REF}"
echo "Reference ready: ${REF}"
