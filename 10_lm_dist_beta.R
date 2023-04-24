
# Load packages
library(cellssm)
library(ggplot2)
library(patchwork)

# Create output directory
out <- "10_lm_dist_beta"
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
mvtime <- as.data.frame(data.table::fread("02_ssm_individual_stepwise=c(99, 90)/csv/ssm_individual_chloroplast_mvtime.csv"))

# robust = F
glist <- lm_dist_beta(cell_list = cell_list, mvtime = mvtime, robust = F,
                      ssm_path = "02_ssm_individual_stepwise=c(99, 90)",
                      ssm_method = "Bayes", res_name = "chloroplast",
                      ex_name = "microbeam", unit1 = "micrometer", unit2 = "min")

g <- (glist[[1]] + labs(tag = "A")) + (glist[[2]] + labs(tag = "B")) + (glist[[3]] + labs(tag = "C")) + 
  (glist[[4]] + labs(tag = "D")) + (glist[[5]] + labs(tag = "E")) + 
  plot_layout(nrow = 2)
suppressWarnings(ggsave(paste0(out, "/lm_dist_beta_robustF.pdf"),
                        g, height = 104, width = 168, units = "mm"))

# robust = T
glist <- lm_dist_beta(cell_list = cell_list, mvtime = mvtime, robust = T,
                      ssm_path = "02_ssm_individual_stepwise=c(99, 90)",
                      ssm_method = "Bayes", res_name = "chloroplast",
                      ex_name = "microbeam", unit1 = "micrometer", unit2 = "min")

g <- (glist[[1]] + labs(tag = "A")) + (glist[[2]] + labs(tag = "B")) + (glist[[3]] + labs(tag = "C")) + 
  (glist[[4]] + labs(tag = "D")) + (glist[[5]] + labs(tag = "E")) + 
  plot_layout(nrow = 2)
suppressWarnings(ggsave(paste0(out, "/lm_dist_beta_robustT.pdf"),
                        g, height = 104, width = 168, units = "mm"))

