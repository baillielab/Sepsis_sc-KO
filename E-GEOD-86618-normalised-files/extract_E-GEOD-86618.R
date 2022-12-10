# 
# Author: Nureen Zaki
# Date: 2022-09-08
# Title: Sepsis-MAIC single cell KO using scTenifoldKnk
#

# Libraries ---------------------------------------------------------------

library(Matrix) # 1.5-1
library(Seurat) # 3.2.3
library(tidyverse)
library(tictoc)

# Directories -------------------------------------------------------------

wd <- 'C:/Users/nuree/Downloads/INTERN_HOME/sc-KO/sepsis/'
scdir <- 'E-GEOD-86618-normalised-files/'
setwd(paste0(wd, scdir))

fle <- 'E-GEOD-86618.aggregated_filtered_normalised_counts'

# Load .mtx file ----------------------------------------------------------

# Make new feature (row) file
row_file <- paste0(fle, '.mtx_rows')

feat_row <- read.csv(row_file, header = FALSE, sep = '\t') %>% dplyr::select(V1)

write.table(feat_row, file = paste0(row_file, '2'), 
            quote = FALSE, row.names = FALSE, col.names = FALSE)

remove(feat_row)

# Filter dataset: get only the normal@WT samples --------------------------

# Get cell IDs of normal/WT samples from exp design filegre
exp_design <- read.csv('ExpDesign-E-GEOD-86618.tsv', sep = '\t') %>%
  filter(Sample.Characteristic.disease. == 'normal')

# Load expression matrix
expression_matrix <- ReadMtx(
  mtx = paste0(sfle, '.mtx'),
  features = paste0(fle, '.mtx_rows2'),
  cells = paste0(fle, '.mtx_cols'),
  feature.column = 1,
  cell.sep = '\n',
  feature.sep = '\n')

# Subset matrix to cells from only 'normal' or 'WT' samples
expression_matrix <- expression_matrix[, colnames(expression_matrix) %in% exp_design$Assay]


# Convert ENSG to Gene Symbol ---------------------------------------------

library(limma) # 3.52.2
library(data.table) # 1.14.2
library(dplyr)
library(stringr)
library(gprofiler2) # 0.2.1
library(biomaRt) # 2.52.0

# Approach 1
mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
genes <- rownames(expression_matrix)
G_list <- getBM(
  filters = "ensembl_gene_id", 
  attributes = c("ensembl_gene_id","hgnc_symbol"),
  values = genes,
  mart = mart)
  ## Does not map all ENSG ids to gene symbol: 26251/26267

# Check which IDs are not mapped to a HGNC symbol by biomaRt using gprofiler2
get_HGNC <- function(x){
  symbol <- gconvert(x, organism = 'hsapiens', target = 'HGNC', filter_na = FALSE) %>% 
    dplyr::select(target) %>% unlist() %>% unname()
  if (is.na(symbol) == TRUE) {symbol <- x}
  return(symbol)
}

not_mapped <- data.frame(ENSG = setdiff(rownames(expression_matrix), G_list$ensembl_gene_id))
not_mapped <- not_mapped %>% mutate(HGNC = lapply(ENSG, get_HGNC))


# Assign feature/column name to expression matrix -------------------------

# Subset the expression matrix by genes that have biomaRt annotation
exp_mat_subset <- expression_matrix[rownames(expression_matrix) %in% G_list$ensembl_gene_id, ]

# Rename rows ENSG -> HGNC
rownames(exp_mat_subset) <- G_list$hgnc_symbol                          

# Save dgCMatrix as RObject
exp_mat_subset <- as.data.frame(exp_mat_subset)
exp_mat_subset <- subset(exp_mat_subset, rowSums(exp_mat_subset)>0)

saveRDS(exp_mat_subset, file = paste0(fle, '_RData.RData'))


write.table(exp_mat_subset, paste0(fle, '_matrix.csv'), quote = FALSE, sep = '\t', row.names = TRUE, col.names = NA)

# QC as seurat object -----------------------------------------------------

min_cells <- .05*dim(exp_mat_subset)[2]
seurat_object <- CreateSeuratObject(counts = exp_mat_subset,
                                    min.cells = min_cells,
                                    min.features = 0)



# seurat_object[['percent.mt']] <- PercentageFeatureSet(seurat_object, pattern = '^MT-')
#   ## No mitochondrial genes
# View(seurat_object@meta.data)
VlnPlot(seurat_object, features = c('nFeature_RNA', 'nCount_RNA'), ncol = 3)
FeatureScatter(seurat_object, feature1 = 'nCount_RNA', feature2 = 'nFeature_RNA') +
  geom_smooth(method = 'lm')

# get_ENSG <- function(x){
#   ensg <- gconvert(x, organism = 'hsapiens', target = 'ENSG', filter_na = FALSE) %>% 
#     dplyr::select(target) %>% unlist() %>% unname()
#   if (is.na(ensg) == TRUE){
#     ensg <- x
#   }
#   return(ensg)
# }

# Dump --------------------------------------------------------------------

# exp_mat <- Matrix::readMM('E-GEOD-86618-quantification-raw-files/E-GEOD-86618.aggregated_filtered_counts.mtx')
# exp_rows <- read.csv('E-GEOD-86618-quantification-raw-files/E-GEOD-86618.aggregated_filtered_counts.mtx_rows',
#                      header = FALSE, sep = '\t') %>% select(V1) %>% unlist()
# exp_cols <- read.csv('E-GEOD-86618-quantification-raw-files/E-GEOD-86618.aggregated_filtered_counts.mtx_cols',
#                      header = FALSE, sep = '\n') %>% unlist()
# 
# dimnames(exp_mat) <- list(exp_rows, exp_cols)

# Approach 2
# library('EnsDb.Hsapiens.v79') # 2.99.0
# geneIDs <- ensembldb::select(EnsDb.Hsapiens.v79, 
#                              keys = exp_rows, 
#                              keytype = "GENEID", 
#                              columns = c("SYMBOL","GENEID")
#                              )
# Does not map all ENSG ids to gene symbol 25296/26267

# Approach 3 // FAILED
# library('org.Hs.eg.db') # 3.15.0
# infofile <- read.csv('../Homo_sapiens.gene_info.txt', sep = '\t')
# find_alias <- function(x){
#   alias <- limma::alias2SymbolUsingNCBI(
#     x, 
#     infofile, 
#     required.columns = c('Symbol', 'Synonyms', 'dbXrefs'))
#   return(alias)
# }
# alias <- lapply(exp_rows, find_alias)

# Approach 4
# library(gprofiler2)
# genes_gprof <- gconvert(query = exp_rows, organism = "hsapiens",
#                         target = 'HGNC', mthreshold = Inf, filter_na = TRUE,
#                         numeric_ns = 'ENSG') %>% dplyr::select(input, target)
# 21298 of 26267 genes mapped
