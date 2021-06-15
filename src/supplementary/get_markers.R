library(Seurat)
library(SeuratData)
library(fcoex)
library(readr)
library(patchwork)
library(ggplot2)
library(dplyr)
library(tidyr)
library(scales)

fc <- read_rds("./data/intermediate_data/fc_pbmc3k.rds")
data("pbmc3k.final")

# Add modules to Seurat Object

module_idents <- data.frame(fc@mod_idents)

colnames(module_idents) <- paste("fcoex", colnames(module_idents), sep="_")

metadata = cbind(pbmc3k.final@meta.data, module_idents)

pbmc3k.final@meta.data <- metadata

markers_m5 <- FindMarkers(pbmc3k.final,
                          ident.1 = "HP",
                          group.by = "fcoex_HLA.DRB1",
                          only.pos = TRUE) %>%
  tibble::rownames_to_column(var = "gene")

markers_m2 <- FindMarkers(pbmc3k.final,
                          ident.1 = "HP",
                          group.by = "fcoex_CD3D",
                          only.pos = TRUE)  %>% 
  tibble::rownames_to_column(var = "gene")

markers_m7 <- FindMarkers(pbmc3k.final,
                          ident.1 = "HP",
                          group.by = "fcoex_CST7",
                          only.pos = TRUE) %>% 
  tibble::rownames_to_column(var = "gene")

markers_m9 <- FindMarkers(pbmc3k.final,
                          ident.1 = "HP",
                          group.by = "fcoex_FCGR3A",
                          only.pos = TRUE) %>% 
  tibble::rownames_to_column(var = "gene")

write_tsv(markers_m5,
          file = "results/tables/markers_hla.tsv" )

write_tsv(markers_m2,
          file = "results/tables/markers_cd3d.tsv" )

write_tsv(markers_m7,
          file = "results/tables/markers_cst7.tsv" )

write_tsv(markers_m9,
          file = "results/tables/markers_fcgr3a.tsv" )


