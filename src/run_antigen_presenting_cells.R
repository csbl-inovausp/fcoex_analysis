library(Seurat)
library(fcoex)
library(readr)
source("src/helpers.R")

filename <-
  "/home/lubianat/Documents/lab_related/fcoex_analysis/data/apc_covid/GSE169346_ValidationSet_RPMI_SeuratObj.Rds"
apc_validation <- readRDS(filename)

# Run fcoex for healthy subset
healthy <-
  subset(x = apc_validation, subset = severity == "Healthy")
rm(apc_validation)
gc()

# Sample for local processing
healthy <-
  healthy[, sample(colnames(healthy), size = 1000, replace = F)]

exprs <-
  as.data.frame(Seurat::GetAssayData(object = healthy, slot = "data"))
target <- healthy@meta.data$celltype

fc <- new_fcoex(data.frame(exprs), target)
fc <- discretize(fc)
fc <-
  find_cbf_modules(
    fc,
    verbose = TRUE,
    n_genes_selected_in_first_step = 100,
    is_parallel = TRUE
  )
fc <- recluster(fc)

write_rds(x = fc,
          file = "./data/intermediate_data/apc_healthy_small.rds",
          compress = c("gz"))


proportions_df = get_proportions_df(seurat_object = healthy,
                   fcoex_object = fc ,
                   ident_column = "celltype",
                   algorithm_name = "fcoex")

ggsave("results/figures/fcoex_apc_healthy.pdf", width = 8, height = 8)
plot_proportions_df(proportions_df, change_level_order_for_plotting = FALSE)
while (!is.null(dev.list()))  dev.off()



# Process fcoex for severe cases
apc_validation <- readRDS(filename)
severe <- subset(x = apc_validation, subset = severity == "Severe")
severe <- severe[, sample(colnames(severe), size = 1000, replace = F)]  

exprs <-
  as.data.frame(Seurat::GetAssayData(object = severe, slot = "data"))
target <- severe@meta.data$celltype
  
fc <- new_fcoex(data.frame(exprs), target)
fc <- discretize(fc)

fc <-
  find_cbf_modules(
    fc,
    verbose = TRUE,
    n_genes_selected_in_first_step = 50,
    is_parallel = TRUE
  )
fc <- recluster(fc)

proportions_df <- get_proportions_df(seurat_object = severe,
                                    fcoex_object = fc ,
                                    ident_column = "celltype",
                                    algorithm_name = "fcoex")

ggsave("results/figures/fcoex_apc_severe.pdf", width = 8, height = 8)
plot_proportions_df(proportions_df, change_level_order_for_plotting = FALSE)
while (!is.null(dev.list()))  dev.off()


ggsave("results/figures/umap_apc_healthy.pdf", width = 8, height = 8)
p1 <- Seurat::DimPlot(healthy, group.by = "celltype") +
  ggtitle("Healthy patients")
plot(p1)
while (!is.null(dev.list()))  dev.off()


ggsave("results/figures/umap_apc_severe.pdf", width = 8, height = 8)
p2 <- Seurat::DimPlot(severe, group.by = "celltype") +
  ggtitle("Severe COVID-19 patients")
plot(p2)
while (!is.null(dev.list()))  dev.off()

write_rds(x = fc,
            file = "./data/intermediate_data/apc_severe_covid.rds",
            compress = c("gz"))
  