library(Seurat)
library(fcoex)


filename <- "/home/lubianat/Documents/lab_related/fcoex_analysis/data/apc_covid/GSE169346_ValidationSet_RPMI_SeuratObj.Rds"
apc_validation <- readRDS(filename)

colnames(apc_validation@meta.data)
Seurat::DimPlot(apc_validation, group.by = "SCT_snn_res.1.2" )

exprs <- as.data.frame(Seurat::GetAssayData(object = apc_validation, slot = "data"))
exprs <- exprs[1:10000, 1:10000]
target < apc_validation@meta.data$SCT_snn_res.1.2[1:10000]


fc <- new_fcoex(data.frame(exprs),target)
fc <- discretize(fc)

fc <- find_cbf_modules(fc)
fc <- recluster(fc)

write_rds(
  x=fc,
  file="./data/intermediate_data/prototype_apc_validation.rds",
  compress = c("gz"))



