# Helper functions
library(tidyr)
library(dplyr)
library(ggplot2)
get_proportions_df <-
  function(seurat_object,
           fcoex_object,
           algorithm_name = "fcoex",
           ident_column = "seurat_annotations") {
    module_idents <- data.frame(fcoex_object@mod_idents)
    colnames(module_idents) <-
      paste(algorithm_name, colnames(module_idents), sep = "_")
    
    # Caution! cbind assumes order is kept the same, what is generally true,
    # but might lead to silent errors if not.
    metadata = cbind(seurat_object@meta.data, module_idents)
    
    metadata <- droplevels(metadata)
    
    
    original_clusters <- data.frame(table(metadata[ident_column]))
    seurat_to_fcoex_proportions = list()
    
    for (identity in levels(metadata[,ident_column])) {

      identity_metadata = metadata[metadata[ident_column] == identity, ]
      
      proportions_for_this_identity = list()
      
      for (col_name in colnames(identity_metadata)) {
        if (grepl(algorithm_name, col_name)) {
          total_cells = nrow(identity_metadata)
          header_positive_cells = sum(identity_metadata[col_name] == "HP")
          proportions_for_this_identity[[col_name]] <-
            header_positive_cells / total_cells
        }
        seurat_to_fcoex_proportions[[identity]] <-
          proportions_for_this_identity
      }
    }
    
    proportions_df <-
      data.frame(do.call(rbind, seurat_to_fcoex_proportions))
    return(proportions_df)
  }

prepare_proportions_df <- function(proportions_df, change_level_order_for_plotting = TRUE) {
  # Cluster to get factor order ----------
  seurat_distances <- dist(proportions_df)
  fcoex_distances <- dist(t(proportions_df))
  seurat_clustering <- hclust(seurat_distances)
  fcoex_clustering <- hclust(fcoex_distances)
  
  # Pivot to longer format -----------
  proportions_df$Seurat <- rownames(proportions_df)
  proportions_df_long <- proportions_df %>%
    pivot_longer(cols = !Seurat,
                 names_to = "fcoex",
                 values_to = "proportion")
  
  proportions_df_long$proportion <-
    as.numeric(proportions_df_long$proportion)
  if (change_level_order_for_plotting){
  # Adjust types and order of levels of columns in dataframe ----------

    proportions_df_long$Seurat <-
      factor(proportions_df_long$Seurat, levels = seurat_clustering$labels[seurat_clustering$order])
    
    proportions_df_long$fcoex <-
      factor(proportions_df_long$fcoex, levels = fcoex_clustering$labels[fcoex_clustering$order])
  }
  proportions_df_long
  
}

plot_proportions_df <- function(proportions_df, change_level_order_for_plotting = TRUE) {
  proportions_df %>%
    prepare_proportions_df(change_level_order_for_plotting) %>%
    ggplot(aes_string(x = "Seurat", y = "fcoex", fill = "proportion")) +
    geom_tile() +
    scale_fill_gradient(low = "white", high = "purple3") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 70, hjust = 1))
  
}