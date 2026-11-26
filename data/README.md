# Input data

Raw inputs are public and are not stored in Git. Download them with:

```bash
bash data/download_data.sh            # everything (~8 GB, SEA-AD needs the AWS CLI)
SKIP_SEAAD=1 bash data/download_data.sh   # GEO datasets only (~48 MB)
```

Expected layout:

```
data/
├── GSE132903/GSE132903_matrix.txt.gz            # GEO series matrix (log2 intensities, 97 AD / 98 ND)
├── GSE104704/GSE104704_RNA-Seq_Table.txt.gz     # GEO supplementary table (size-factor-normalised counts)
├── GSE104704/GSE104704_series_matrix.txt.gz     # sample annotations
└── SEA-AD-MTG/donors_objects/*.h5ad             # 20 SEA-AD MTG donors (10 dementia / 10 no dementia)
```

| Dataset | Source | Notes |
|---|---|---|
| GSE132903 | [NCBI GEO](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE132903) | Illumina HumanHT-12 v4, middle temporal gyrus |
| GSE104704 | [NCBI GEO](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE104704) | RNA-seq, lateral temporal lobe; raw integer counts are recovered exactly from the normalised table in script 02 |
| SEA-AD | [Seattle Alzheimer’s Disease Brain Cell Atlas](https://portal.brain-map.org/explore/seattle-alzheimers-disease) | `s3://sea-ad-single-cell-profiling/MTG/RNAseq/donors_objects/`; raw counts in `layers["UMIs"]` |

The STRING v12.5 and TRRUST v2 files used by the network analysis are small and are included in `network/data/`.
