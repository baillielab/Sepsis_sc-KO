# 
# Author: Nureen Zaki
# Date: 2022-09-08
# Title: Sepsis-MAIC single cell KO using scTenifoldKnk, data from Broad Institute Single Cell Portal
#

# Libraries ---------------------------------------------------------------

library(Matrix) # 1.5-1
library(Seurat) # 3.2.3
library(tidyverse)

# Directories -------------------------------------------------------------

wd <- 'C:/Users/nuree/Downloads/INTERN_HOME/sc-KO/sepsis/'
setwd(wd)

scdir <- 'SCP1219/'
fle <- 'lung'

# Load .mtx file ----------------------------------------------------------

# Loading example single-cell data
# scRNAseq <- system.file("single-cell/example.csv", package = "scTenifoldKnk")
# scRNAseq <- read.csv(scRNAseq, row.names = 1)

# Make new feature (row) file
# row_file <- paste0(scdir, fle, '_genes.tsv')
# feat_row <- read.csv(row_file, header = FALSE, sep = '\t') 
# write.table(feat_row, file = paste0(row_file, '2'), 
#             quote = FALSE, row.names = FALSE, col.names = FALSE)
# remove(feat_row)

# Filter dataset: get only the normal@WT samples --------------------------

# Get cell IDs of normal/WT samples from exp design filegre
metadata <- list.files(path = paste0('./',scdir), pattern='.*(meta|design).*\\.(csv|tsv|txt)$')
exp_design <- read.csv(metadata, sep = '\t') %>%
  filter(organ__ontology_label == 'lung') %>%
  filter(disease__ontology_label == 'normal')

# Load expression matrix
mtx <- list.files(path = paste0('./',scdir), pattern='\\.mtx$')
features <- list.files(path = paste0('./',scdir), pattern='.*(gene|row|feature).*\\.(csv|tsv)$')
cells <- list.files(path = paste0('./',scdir), pattern='.*(cell|barcode|column).*\\.(csv|tsv)$')

expression_matrix <- ReadMtx(
  mtx = mtx,
  features = features,
  cells = cells,
  feature.column = 1,
  cell.sep = '\n',
  feature.sep = '\n')

# Subset matrix to cells from only 'normal' or 'WT' samples
expression_matrix <- expression_matrix[, colnames(expression_matrix) %in% exp_design$NAME]


# Convert ENSG to Gene Symbol ---------------------------------------------

library(data.table) # 1.14.2
library(dplyr)
library(stringr)

# Approach 1
if (startsWith(rownames(expression_matrix)[1], 'ENSG')){
  library(biomaRt) # 2.52.0
  library(limma) # 3.52.2
  library(gprofiler2) # 0.2.1
  mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
  genes <- rownames(expression_matrix)
  G_list <- getBM(
    filters = "ensembl_gene_id",
    attributes = c("ensembl_gene_id","hgnc_symbol"),
    values = genes,
    mart = mart)
  
  # Check which IDs are not mapped to a HGNC symbol by biomaRt using gprofiler2
  get_HGNC <- function(x){
    symbol <- gconvert(x, organism = 'hsapiens', target = 'HGNC', 
                       filter_na = FALSE) %>% 
      dplyr::select(target) %>% unlist() %>% unname()
    if (is.na(symbol) == TRUE) {symbol <- x}
    return(symbol)
  }
  
  not_mapped <- data.frame(ENSG = setdiff(rownames(expression_matrix), G_list$ensembl_gene_id))
  not_mapped <- not_mapped %>% mutate(HGNC = lapply(ENSG, get_HGNC))
  
  # Subset the expression matrix by genes that have biomaRt annotation
  # exp_mat_subset <- expression_matrix[rownames(expression_matrix) %in% G_list$ensembl_gene_id, ]
  
  # Rename rows ENSG -> HGNC
  # rownames(exp_mat_subset) <- G_list$hgnc_symbol   
  } else {
  print('IDs are not ENSG. No IDs converted.')
  }


# QC as Seurat object -----------------------------------------------------

seurat_object <- CreateSeuratObject(counts = expression_matrix,
                                    min.cells = 50,
                                    min.features = 2500)

seurat_object[['percent.mt']] <- PercentageFeatureSet(seurat_object, pattern = '^MT-')
## No mitochondrial genes
View(seurat_object@meta.data)
VlnPlot(seurat_object, features = c('nFeature_RNA', 'nCount_RNA', 'percent.mt'), ncol = 3)
FeatureScatter(seurat_object, feature1 = 'nCount_RNA', feature2 = 'nFeature_RNA') +
  geom_smooth(method = 'lm')

# Save dgCMatrix as RObject
exp_mat_subset <- as.data.frame(seurat_object@assays$RNA@data)
saveRDS(exp_mat_subset, file = paste0(scdir, fle, '_RData.RData'))
