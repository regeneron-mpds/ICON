# ICON
ICON is an R package for identifing reliable TCR-pMHC interactions from multi-omics multiplexed multimer binding single cell sequencing data (Zhang et al. 2020).

## Install and load ICON package:
```
library('devtools')
install_github("regeneron-mpds/ICON") # R version >= 3.6.1
library('ICON')
```
## Steps to run ICON:
Download the [Demo Data](https://github.com/regeneron-mpds/ICON/tree/master/tests/testdata/) in this repo for an example run.

**Step 1**: filter out low-quality cells based on scRNA data
```
rna.raw <- readRDS("path_to_demo_data/rna_raw.rds") # input single cell RNA-seq data
adt.raw <- readRDS("path_to_demo_data/adt_raw.rds") # input single cell dextramer data
pQC.data <- QC(rna.raw)
good.cells <- selectGoodCells(pQC.data, adt.raw)
test.data <- good.cells[15:58,] # extract test dextramer data
nc.data <- good.cells[59:64,] # extract negative control dextramer data
```
**Step 2**: select T cells with paired alpha and beta chains
```
tcr.all <- readRDS("path_to_demo_data/filtered_tcr_annotation.rds")
tcr.pair <- pair_TCR(tcr.all)
```
**Step 3**: identify pMHC binding T cells
```
bg <- estBgNoise(test.data, nc.data) # inspect background noise to choose cutoff
binders <- Identify_binders(test.data, tcr.pair, out_dir = path_to_output_files)
```
