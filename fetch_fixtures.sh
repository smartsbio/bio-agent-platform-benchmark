#!/usr/bin/env bash
# Fetch the fixtures that are too large to keep in git, and regenerate the
# uncompressed read files. Safe to re-run; it skips whatever is already present.
set -euo pipefail
cd "$(dirname "$0")/fixtures"

if [ ! -f pbmc3k_raw.h5ad ]; then
  echo "fetching pbmc3k_raw.h5ad (~6 MB) ..."
  curl -fsSL -o pbmc3k_raw.h5ad \
    "https://datasets.cellxgene.cziscience.com/pbmc3k.h5ad" \
    || echo "FAILED — download 10x Genomics PBMC 3k manually and save it here as pbmc3k_raw.h5ad"
else
  echo "pbmc3k_raw.h5ad already present"
fi

for n in 1 2; do
  if [ ! -f "reads_R${n}.fastq" ] && [ -f "reads_R${n}.fastq.gz" ]; then
    echo "decompressing reads_R${n}.fastq.gz ..."
    gunzip -kf "reads_R${n}.fastq.gz"
  fi
done

echo "fixtures ready: $(ls -1 | wc -l | tr -d ' ') files"
