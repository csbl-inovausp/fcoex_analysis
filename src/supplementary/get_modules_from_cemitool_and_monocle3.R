library(CEMiTool)
library(SeuratData)
library(monocle3)
library(dplyr) # imported for some downstream data manipulation
library(fcoex)
library(Seurat)
data("pbmc3k.final")

exprs <- GetAssayData(object = pbmc3k.final, slot = "data")
meta <- data.frame(cell_id = Idents(pbmc3k.final))

annot <-
  data.frame(SampleName = colnames(exprs), Class = meta$cell_id)

# Get modules via CEMiTool ----------

gmt_fname <-
  system.file("extdata", "pathways.gmt", package = "CEMiTool")
gmt_in <- read_gmt(gmt_fname)
int_df <-
  read.delim(system.file("extdata", "interactions.tsv", package = "CEMiTool"))
cem <- cemitool(
  data.frame(exprs),
  annot,
  gmt = gmt_in,
  interactions = int_df,
  verbose = TRUE,
  plot = TRUE
)

gene_module_df_clean <- cem@module
gene_module_df_clean$modules <-
  as.factor(gene_module_df_clean$modules)

cemitool_modules <- list()
for (module_name in levels(gene_module_df_clean$modules)) {
  module_genes <-
    gene_module_df_clean[gene_module_df_clean$modules == module_name,]$genes
  cemitool_modules[[module_name]] <- module_genes
}
saveRDS(cemitool_modules, file = "data/intermediate_data/pbmc3k_modules_cemitool.rds")

rm(list = ls())

# Get modules via monocle3 ----------

data("pbmc3k.final")

exprs <- GetAssayData(object = pbmc3k.final, slot = "data")
meta <- data.frame(cell_id = Idents(pbmc3k.final))

annot <-
  data.frame(SampleName = colnames(exprs), Class = meta$cell_id)


cds <- new_cell_data_set(exprs, meta)
cds <- preprocess_cds(cds, num_dim = 100)
cds <- reduce_dimension(cds)
cds <- cluster_cells(cds)
pr_graph_test_res <- graph_test(cds, neighbor_graph = "knn", cores = 8)
pr_deg_ids <-
  row.names(subset(pr_graph_test_res, morans_I > 0.01 &
                     q_value < 0.05))
gene_module_df <-
  find_gene_modules(cds[pr_deg_ids], resolution = 1e-3)
gene_module_df_clean <- gene_module_df[, c("id", "module")]

monocle3_modules <- list()
for (i in levels(gene_module_df_clean$module)) {
  module_genes <-
    gene_module_df_clean[gene_module_df_clean$module == i,]$id
  module_name = paste0("M", i)
  monocle3_modules[[module_name]] <- module_genes
}

saveRDS(monocle3_modules, file = "data/intermediate_data/pbmc3k_modules_monocle3.rds")
