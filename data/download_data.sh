#!/usr/bin/env bash
# Download the public input data into data/ (run from the repository root).
#   GSE132903  microarray series matrix                     ~42 MB
#   GSE104704  RNA-seq table + series matrix                ~5.5 MB
#   SEA-AD     20 donor h5ad files (MTG snRNA-seq)          ~8 GB  (needs the AWS CLI; no account required)
# SHA-256 checksums of the GEO files used in this project are verified after download.
set -euo pipefail
cd "$(dirname "$0")"

GEO=https://ftp.ncbi.nlm.nih.gov/geo/series
mkdir -p GSE132903 GSE104704 SEA-AD-MTG/donors_objects
curl -fL -o GSE132903/GSE132903_matrix.txt.gz          "$GEO/GSE132nnn/GSE132903/matrix/GSE132903_series_matrix.txt.gz"
curl -fL -o GSE104704/GSE104704_RNA-Seq_Table.txt.gz   "$GEO/GSE104nnn/GSE104704/suppl/GSE104704_RNA-Seq_Table.txt.gz"
curl -fL -o GSE104704/GSE104704_series_matrix.txt.gz   "$GEO/GSE104nnn/GSE104704/matrix/GSE104704_series_matrix.txt.gz"

shasum -a 256 -c <<SUMS
11286533c0730ccecb97ae3cc773cc51f87b53f631ea5e3c4893a1fb6ac45290  GSE132903/GSE132903_matrix.txt.gz
457cf589f67a5aa90e3efde848792c176b97ce69325e4d85fd0e45dd0fd760c8  GSE104704/GSE104704_RNA-Seq_Table.txt.gz
339148fc8a8f77cf02b0d4ca66b898a576b2c9e8cce2b45f81e124863fec3cdf  GSE104704/GSE104704_series_matrix.txt.gz
SUMS

if [[ "${SKIP_SEAAD:-0}" == "1" ]]; then echo "SKIP_SEAAD=1: skipping SEA-AD download"; exit 0; fi
command -v aws >/dev/null || { echo "Install the AWS CLI to download SEA-AD (or set SKIP_SEAAD=1)"; exit 1; }
DONORS="H20.33.004 H20.33.018 H20.33.020 H20.33.028 H20.33.033 H20.33.037 H20.33.038 H20.33.041 H21.33.027 H21.33.029
        H19.33.004 H20.33.001 H20.33.014 H20.33.019 H20.33.025 H20.33.034 H20.33.035 H21.33.003 H21.33.011 H21.33.025"
for d in $DONORS; do
  aws s3 cp --no-sign-request \
    "s3://sea-ad-single-cell-profiling/MTG/RNAseq/donors_objects/${d}_SEAAD_MTG_RNAseq_final-nuclei.2026-06-22.h5ad" \
    SEA-AD-MTG/donors_objects/
done
echo "All inputs downloaded."
