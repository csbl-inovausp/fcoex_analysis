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

coex_module_plots <- show_net(fc)
modules_to_filter <-module_genes(fc)

# Plot heatmaps for each of the selected modules -------

heat_1 <- DoHeatmap(object = pbmc3k.final, 
                    slot = "data",
                    size = 5,
                    raster = FALSE,
                    features = c(modules_to_filter[["HLA-DRB1"]])) +
  scale_fill_gradient(low = "white", high = "black")  +
  theme(legend.position = "none",
        axis.line.y = element_line()) +
  ylab("M5 (HLA-DRB1)")


heat_2 <- DoHeatmap(object = pbmc3k.final, 
                    label = FALSE,
                    slot = "data",
                    size = 5,
                    raster = FALSE,
                    features = c(modules_to_filter[["CD3D"]])) +
  scale_fill_gradient(low = "white", high = "black")  +
  theme(legend.position = "none",
        axis.line.y = element_line()) +
  ylab("M2 (CD3D)")


heat_3 <- DoHeatmap(object = pbmc3k.final, 
                    label = FALSE,
                    slot = "data",
                    size = 5,
                    raster = FALSE,
                    features = c(modules_to_filter[["CST7"]])) +
  scale_fill_gradient(low = "white", high = "black")  +
  theme(legend.position = "none",
        axis.line.y = element_line()) +
  ylab("M8 (CST7)")

p_heat <- (heat_1 / heat_2 / heat_3)


# Plot the reclusterings for each of the selected modules -------

# Add modules to Seurat Object

module_idents <- data.frame(fc@mod_idents)

colnames(module_idents) <- paste("fcoex", colnames(module_idents), sep="_")

metadata = cbind(pbmc3k.final@meta.data, module_idents)

pbmc3k.final@meta.data <- metadata



recluster_theme <- function(){
  theme(#axis.title=element_blank(),
    #axis.line =element_blank(),
    axis.text=element_blank(),
    axis.ticks=element_blank(),
    legend.position = "none",
    legend.direction = "vertical",
    plot.title = element_text(hjust = 0.5, face = "plain", size = 20))
}

p1 <- DimPlot(pbmc3k.final) +
  theme(# axis.title=element_blank(),
    # axis.line =element_blank(),
    axis.text=element_blank(),
    axis.ticks=element_blank(),
    plot.title = element_text(hjust = 0.5, face = "plain", size = 25))+
  ggtitle("pbmc3k")
p2 <- DimPlot(pbmc3k.final, group.by = "fcoex_HLA.DRB1", cols = c("gray70", "black")) + 
  recluster_theme() +
  ggtitle("C5 (HLA-DRB1)")
p3 <- DimPlot(pbmc3k.final, group.by = "fcoex_CD3D", cols = c("gray70", "black")) +
  recluster_theme() +
  ggtitle("C2 (CD3D)")
p4 <- DimPlot(pbmc3k.final, group.by = "fcoex_CST7", cols = c("gray70", "black")) +
  recluster_theme() +
  ggtitle("C8 (CST7)")

layout <- "
AABBBCCDDEE
AABBBCCDDEE
"

ggsave("results/figures/fcoex_pipeline_pbmc3k.pdf", width = 20, height = 7)
p1 + p_heat + p2 + p3 +p4 + plot_layout(design = layout) 
dev.off()

original_clusters <- data.frame(table(Idents(pbmc3k.final)))

metadata <- pbmc3k.final@meta.data

seurat_to_fcoex_proportions = list()
for (identity in levels(Idents(pbmc3k.final))){
  
  identity_metadata = metadata[metadata["seurat_annotations"] == identity,]
  
  proportions_for_this_identity = list()
  
  for (col_name in colnames(identity_metadata)){
    if (grepl("fcoex", col_name)){
      
      total_cells = nrow(identity_metadata)
      header_positive_cells = sum(identity_metadata[col_name] == "HP")
      proportions_for_this_identity[[col_name]] <- header_positive_cells/total_cells 
    }
    seurat_to_fcoex_proportions[[identity]] <- proportions_for_this_identity
  }
}

proportions_df <- data.frame(do.call(rbind, seurat_to_fcoex_proportions))

# Plot the proportions in relation to the Seurat object ----------


prepare_proportions_df <-function(proportion_df){
  # Cluster to get factor order ----------
  seurat_distances <- dist(proportions_df)
  fcoex_distances <- dist(t(proportions_df))
  seurat_clustering <- hclust(seurat_distances)
  fcoex_clustering <- hclust(fcoex_distances)
  
  # Pivot to longer format -----------
  proportions_df$Seurat <- rownames(proportions_df)
  proportions_df_long <- proportions_df %>%
    pivot_longer(cols = !Seurat,names_to="fcoex", values_to="proportion") 
  
  # Adjust types and order of levels of columns in dataframe ----------
  proportions_df_long$proportion <- as.numeric(proportions_df_long$proportion)
  
  proportions_df_long$Seurat <- factor(proportions_df_long$Seurat, levels = seurat_clustering$labels[seurat_clustering$order])
  
  proportions_df_long$fcoex <- factor(proportions_df_long$fcoex, levels = fcoex_clustering$labels[fcoex_clustering$order])
  
  proportions_df_long
  
}

ggsave("results/figures/fcoex_reclusterings_pbmc3k.pdf", width = 4, height= 4)
proportions_df %>%
  prepare_proportions_df() %>%
  ggplot(aes_string(x="Seurat", y="fcoex", fill="proportion")) +
  geom_tile() +
  scale_fill_gradient(low = "white", high = "purple3") +
  theme_bw()+
  theme(axis.text.x = element_text(angle = 70, hjust = 1))
dev.off()


