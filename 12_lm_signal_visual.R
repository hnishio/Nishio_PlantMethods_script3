
# Load packages
library(cellssm)
library(ggplot2)
library(patchwork)

# Create output directory
out <- "12_lm_signal_visual"
if(file.exists(out)==F){
  dir.create(out, recursive=T)
}

# Load data
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
visual <- read.csv("data/visual_judgement.csv")
mvtime <- as.data.frame(data.table::fread("04_nomodel/nomodel_chloroplast_mvtime.csv"))
mvtime <- mvtime[,-4]

# robust = F
glist <- lm_signal(cell_list = cell_list, mvtime = mvtime, robust = F,
                   ex_name = "microbeam", unit1 = "micrometer", unit2 = "min")

g <- (glist[[1]] + labs(tag = "A")) + (glist[[2]] + labs(tag = "B")) + (glist[[3]] + labs(tag = "C")) + 
  (glist[[4]] + labs(tag = "D")) + (glist[[5]] + labs(tag = "E")) + (glist[[6]] + labs(tag = "F")) + 
  (glist[[7]] + labs(tag = "G")) + (glist[[8]] + labs(tag = "H")) + (glist[[9]] + labs(tag = "I")) + 
  (glist[[10]] + labs(tag = "J")) + (glist[[11]] + labs(tag = "K")) +
  plot_layout(ncol = 3)
suppressWarnings(ggsave(paste0(out, "/lm_signal_visual_robustF.pdf"),
                        g, height = 220, width = 50*3, units = "mm"))

# robust = T
glist <- lm_signal(cell_list = cell_list, mvtime = mvtime, robust = T,
                   ex_name = "microbeam", unit1 = "micrometer", unit2 = "min")

g <- (glist[[1]] + labs(tag = "A")) + (glist[[2]] + labs(tag = "B")) + (glist[[3]] + labs(tag = "C")) + 
  (glist[[4]] + labs(tag = "D")) + (glist[[5]] + labs(tag = "E")) + (glist[[6]] + labs(tag = "F")) + 
  (glist[[7]] + labs(tag = "G")) + (glist[[8]] + labs(tag = "H")) + (glist[[9]] + labs(tag = "I")) + 
  (glist[[10]] + labs(tag = "J")) + (glist[[11]] + labs(tag = "K")) +
  plot_layout(ncol = 3)
suppressWarnings(ggsave(paste0(out, "/lm_signal_visual_robustT.pdf"),
                        g, height = 220, width = 50*3, units = "mm"))

