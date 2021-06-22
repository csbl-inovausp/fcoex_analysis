library(fcoex)
library(patchwork)
library(ggplot2)
library(ggraph)
library(ggrepel)

fc <- readRDS("./data/intermediate_data/fc_pbmc3k.rds")
fcoex_modules <- fc@mod_idents

# Plot ggraph visualization of fcoex networks ---------

fc <- get_nets(fc, n = 8, min_elements = 1)
list_of_network_plots  = list()

temp <- show_net(fc)
for (i in 1:length(temp)) {
  print(i)
  list_of_network_plots[[i]] =  temp[[i]]
}

ggsave("results/supplementary/fcoex_nets.pdf",
       width = 14,
       height = 14)
wrap_plots(list_of_network_plots) + plot_annotation(tag_levels = 'A')
dev.off()

# Plot Over Representation Analysis of fcoex networks ---------

gmt_object <-
  pathwayPCA::read_gmt("data/ReactomePathways.gmt",  description = TRUE)

for (i in gmt_object$pathways) {
  print(i)
}

fc <- mod_ora(fc, gmt_object)
fc <- plot_ora(fc, n = 4, pv_cut = 0.05)

list_of_ora_plots  = list()
for (i in names(fcoex_modules)) {
  list_of_ora_plots[[i]] = fc@barplot_ora[[i]]$pl
}

ggsave("results/supplementary/fcoex_ora.pdf",
       width = 14,
       height = 9)
wrap_plots(list_of_ora_plots) + plot_annotation(tag_levels = 'A')
dev.off()
