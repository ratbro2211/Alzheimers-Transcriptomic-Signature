#!/usr/bin/env bash
# Reproducibility check: rerun the steps that need only the GEO data (~48 MB, no SEA-AD)
# and confirm their result tables are identical to the committed ones.
#   bash data/download_data.sh  (or SKIP_SEAAD=1 bash data/download_data.sh)  first
#   bash tests/reproduce_check.sh
# Exit code 0 = all committed results reproduced exactly.
set -uo pipefail
cd "$(dirname "$0")/.."
PY="${PYTHON:-python3}"
fail=0
step() { echo "== $1"; shift; "$@" > /dev/null 2>&1 || { echo "   FAILED to run: $*"; fail=1; }; }

step "01 microarray DE"          Rscript analysis/dataset1_microarray/scripts/01_limma_analysis.R
step "01b microarray figures"    Rscript analysis/dataset1_microarray/scripts/01b_figures.R
step "02 RNA-seq DE"             Rscript analysis/dataset2_rnaseq/scripts/02_deseq2_analysis.R
step "02b RNA-seq figures"       Rscript analysis/dataset2_rnaseq/scripts/02b_figures.R
step "04 comparison"             Rscript analysis/comparison/scripts/04_comparison.R
step "05a validate microarray"   Rscript analysis/validation/scripts/05a_validate_microarray.R
step "05b validate RNA-seq"      Rscript analysis/validation/scripts/05b_validate_rnaseq.R
step "p2_01 STRING network"      bash -c "cd network && $PY scripts/p2_01_string_network.py"

echo "== comparing result tables with the committed versions"
CHECK=(
  analysis/dataset1_microarray/results/microarray_DEGs.csv
  analysis/dataset2_rnaseq/results/rnaseq_DEGs.csv
  analysis/comparison/results/pairwise_correlations.csv
  analysis/comparison/results/panel_33gene_summary.csv
  analysis/validation/results/checks_microarray.csv
  analysis/validation/results/checks_rnaseq.csv
  analysis/validation/results/loo_rnaseq.csv
  network/results/hub_genes_string.csv
  network/results/network_null_model.csv
  network/results/louvain_communities.csv
)
for f in "${CHECK[@]}"; do
  if git diff --quiet -- "$f"; then echo "   identical  $f"; else echo "   DIFFERENT  $f"; fail=1; fi
done
[[ $fail -eq 0 ]] && echo "PASS: all committed results reproduced" || echo "FAIL: see above"
exit $fail
