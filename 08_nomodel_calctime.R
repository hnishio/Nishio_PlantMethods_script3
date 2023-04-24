
# Load packages
library(cellssm)
library(tictoc)

# Create output directory
out <- "08_nomodel_calctime"
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
cell_list <- list(cell2)
df_cor <- as.data.frame(data.table::fread(file = paste0("04_nomodel", "/grid_cor_mae_rmse.csv")))
df_cor <- df_cor[df_cor$no_data>100,]
df_cor_best3 <- df_cor[which.min(df_cor$rmse),]
df_cor_best <- df_cor_best3

# Predict movement
tic()
nomodel(cell_list = cell_list, out = out,
        res_name = "chloroplast", ex_name = "microbeam", 
        graph = F, unit1 = "micrometer", unit2 = "min",
        consecutive = df_cor_best$consecutive, 
        period = df_cor_best$period, fold = df_cor_best$fold)
tictoc_time <- toc()
calc_time <- data.frame(
  calc_time = as.numeric(tictoc_time$toc - tictoc_time$tic))
data.table::fwrite(calc_time, file = paste0(out, "/nomodel_calctime_cell2.csv"))

