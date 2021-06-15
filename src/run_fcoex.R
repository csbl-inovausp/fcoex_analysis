library(Seurat)
library(SeuratData)
library(fcoex)
library(readr)

data("pbmc3k.final")

exprs <- as.data.frame(GetAssayData(object = pbmc3k.final, slot = "data"))
target <- Idents(pbmc3k.final)


fc <- new_fcoex(data.frame(exprs),target)
fc <- discretize(fc)

fc <- find_cbf_modules(fc)
fc <- recluster(fc)

write_rds(
  x=fc,
  file="./data/intermediate_data/fc_pbmc3k.rds",
  compress = c("gz"))





