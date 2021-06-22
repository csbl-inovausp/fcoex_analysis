# Compare the different module-building tools for single-cell RNA-seq

library(fcoex)
library(Seurat)
library(SeuratData)
library(gridExtra)
library(scales)
library(patchwork)
source("src/helpers.R")

monocle3_modules <-
  readRDS(file = "data/intermediate_data/pbmc3k_modules_monocle3.rds")
cemitool_modules <-
  readRDS(file = "data/intermediate_data/pbmc3k_modules_cemitool.rds")
fc <- readRDS("./data/intermediate_data/fc_pbmc3k.rds")
fcoex_modules <- fc@module_list

data("pbmc3k.final")
# Add CEMiTool modules to fcoex object and recluster ------------

## Only using modules with 2 or more genes -----------
cemitool_modules <-
  cemitool_modules[lapply(cemitool_modules, length) > 2]

fc_cemitool <- fc
fc_cemitool@module_list <- cemitool_modules

##  Assign a random header -----------
names(fc_cemitool@module_list) <- lapply(cemitool_modules, head,1)

fc_cemitool <- recluster(fc_cemitool)

##  Substitute names back -----------
names(fc_cemitool@module_list) <- lapply(cemitool_modules, head,1)
  lapply("CEMiTool_M", paste0, 1:length(cemitool_modules))[[1]]

proportions_df_cemitool <-
  get_proportions_df(pbmc3k.final, fc_cemitool, algorithm_name = "CEMiTool")


ggsave(
  "results/supplementary/cemitool_reclusterings_pbmc3k.pdf",
  width = 4,
  height = 4
)
plot_proportions_df(proportions_df_cemitool) + ylab("CEMiTool")

dev.off()

# Add monocle3 modules to fcoex object and recluster ------------


fc_monocle <- fc
fc_monocle@module_list <- monocle3_modules

##  Assign a random header -----------
names(fc_monocle@module_list) <- lapply(monocle3_modules, head,1)

fc_monocle <- recluster(fc_monocle)

##  Substitute names back -----------
names(fc_monocle@module_list) <-
  lapply("monocle3_M", paste0, 1:length(monocle3_modules))[[1]]


proportions_df_monocle <-
  get_proportions_df(pbmc3k.final, fc_monocle, algorithm_name = "monocle3")

ggsave(
  "results/supplementary/monocle3_reclusterings_pbmc3k.pdf",
  width = 4,
  height = 4
)


# Plot comparison of proportions in relation to original labels --------
plot_proportions_df(proportions_df_monocle) + ylab("monocle3")
dev.off()


proportions_df_fcoex <-
  get_proportions_df(pbmc3k.final, fc, algorithm_name = "fcoex")


p1 <- plot_proportions_df(proportions_df_cemitool) + ylab("CEMiTool") + theme(axis.title.y = element_blank())
p2 <- plot_proportions_df(proportions_df_monocle) + ylab("monocle3") + theme(axis.title.y  = element_blank())
p3 <- plot_proportions_df(proportions_df_fcoex) + ylab("fcoex") + theme(axis.title.y  = element_blank())


ggsave(
  "results/supplementary/comparative_reclusterings_pbmc3k.pdf",
  width = 12,
  height = 4
)
p1 + p2 + p3 + plot_annotation(tag_levels = 'A')
dev.off()


# Plot module sizes --------

module_sizes <- c(unlist(lapply(cemitool_modules, length)),
  unlist(lapply(monocle3_modules,length)),
  unlist(lapply(fcoex_modules,length))
)


module_origin <- factor(c(rep("CEMiTool", length(cemitool_modules)),
                   rep("Monocle3", length(monocle3_modules)),
                   rep("Fcoex", length(fcoex_modules))), levels=c("CEMiTool", "Monocle3", "Fcoex")
)

p4 <- data.frame(n=module_sizes, x=module_origin) %>% 
  ggplot(aes(y=n, x=x, fill=x)) + 
  geom_boxplot() +
  ylab("Number of genes per module") +
  xlab("Algorithm") +
  scale_fill_manual(values= c("beige", "pink", "turquoise")) +
  theme_minimal() +
  theme(legend.position = "none") +
  theme(axis.text.x = element_text(angle = 70, hjust = 1))


ggsave( 
  "results/supplementary/pbmc3k_genes_per_module_count.pdf",
  width = 12,
  height = 8
)
p4  +  theme(legend.position = "left")
dev.off()


layout <- '
ABC
D##
'
wrap_plots(p1, p2, p3, p4, design = layout) + plot_annotation(tag_levels = 'A')
dev.off()



ggsave( 
  "results/supplementary/comparative_reclusterings_pbmc3k_with_count.pdf",
  width = 12,
  height = 8
)

layout <- '
ABC
D##
'
wrap_plots(p1, p2, p3, p4, design = layout) + plot_annotation(tag_levels = 'A')
dev.off()



