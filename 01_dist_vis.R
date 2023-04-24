
# Load packages
library(cellssm)
library(ggplot2)
library(patchwork)

# Create an output directory
out <- "01_dist_vis"
if(file.exists(out)==F){
  dir.create(out, recursive=T)
}

# Load data of chloroplasts
cell1 <- read.csv("data/WT_cell1_v4.csv")
cell2 <- read.csv("data/cell2.csv")
cell3 <- read.csv("data/cell3.csv")
cell4 <- read.csv("data/cell4.csv")
cell5 <- read.csv("data/WT_cell5_v3.csv")
cell6 <- read.csv("data/WT_cell6_v3.csv")
cell7 <- read.csv("data/WT_cell7_v3.csv")
cell8 <- read.csv("data/WT_cell8_v3.csv")
cell9 <- read.csv("data/WT_cell9_v3.csv")
cell_list <- list(cell1, cell2, cell3, cell4, cell5, cell6, cell7, cell8, cell9)

# Plotting
glist <- dist_vis(cell_list = cell_list,
                  res_name = "chloroplast", ex_name = "microbeam", 
                  unit1 = "micrometer", unit2 = "min")

void <- ggplot() + theme_void()
g <- (void + labs(tag = "A")) +
  (glist[[1]] + labs(tag = "B")) + (glist[[2]] + labs(tag = "C")) + (glist[[3]] + labs(tag = "D")) + 
  (glist[[4]] + labs(tag = "E")) + (glist[[5]] + labs(tag = "F")) + (glist[[6]] + labs(tag = "G")) + 
  (glist[[7]] + labs(tag = "H")) + (glist[[8]] + labs(tag = "I")) + (glist[[9]] + labs(tag = "J")) +
  plot_layout(ncol = 2) &
  theme(plot.tag = element_text(size = 12, face = "bold"))
ggsave(paste0(out, "/dist_vis.pdf"),
       g, height = 220, width = 160, units = "mm")

