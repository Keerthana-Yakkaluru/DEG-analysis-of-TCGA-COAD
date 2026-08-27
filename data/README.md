# Data

Raw RNA-seq data is **not tracked** in this repository due to file size (522 samples of STAR-Counts gene expression data from GDC).

## How to regenerate

Run the data acquisition step from `scripts/tcga_coad_deseq2_analysis.R`, which uses `TCGAbiolinks` to query and download the data directly from the GDC (Genomic Data Commons) database:

```r
library(TCGAbiolinks)

query <- GDCquery(
  project = "TCGA-COAD",
  data.category = "Transcriptome Profiling",
  data.type = "Gene Expression Quantification",
  workflow.type = "STAR - Counts",
  sample.type = c("Primary Tumor", "Solid Tissue Normal")
)

GDCdownload(query)
data <- GDCprepare(query)
```

This will download data into a local `GDCdata/` folder (excluded via `.gitignore`) and assemble it into a `RangedSummarizedExperiment` object (60,660 genes × 522 samples: 481 Primary Tumor, 41 Solid Tissue Normal).

**Note:** this step can take a significant amount of time and disk space (several GB) on first run, depending on connection speed.
