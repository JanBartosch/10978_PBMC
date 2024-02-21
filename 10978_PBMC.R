library(dplyr)
library(Seurat)
library(patchwork)

# 1. Open Data and create Seurat Object
# Download data, decompress with Zip7 and store in your desired folder

# Load the PBMC dataset
pbmc.data <- ReadMtx(mtx = "C:/Users/Admin/Desktop/R-script/10978_PBMC/Raw_Data/matrix.mtx",
                     features = "C:/Users/Admin/Desktop/R-script/10978_PBMC/Raw_Data/features.tsv",
                     cells =  "C:/Users/Admin/Desktop/R-script/10978_PBMC/Raw_Data/barcodes.tsv")

# To check, if the cell count matrix is correct, look at the first 10 columns and rows
# columns = barcodes ; rows = features/genes
pbmc.data[1:10,1:10]

# Initialize the Seurat object with the raw (non-normalized data).
pbmc <- CreateSeuratObject(counts = pbmc.data, project = "pbmc3k", min.cells = 3, min.features = 200)
pbmc
str(pbmc)
