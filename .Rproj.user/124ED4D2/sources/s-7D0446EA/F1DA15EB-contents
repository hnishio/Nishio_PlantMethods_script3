
# Load packages
library(cellssm)
library(ggplot2)

# Create output directory
out <- "04_nomodel"
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



### Predict movement (grid search)
grid_par <- expand.grid(consecutive = seq(1, 5, 1),
                        period = seq(1, 15, 1), 
                        fold = seq(1, 4, 0.1))
vec_cor <- NULL
vec_mae <- NULL
vec_rmse <- NULL
no_perfect <- NULL
no_data <- NULL
no_data_full <- NULL
for(i in 1:nrow(grid_par)){
  nomodel(cell_list = cell_list, visual = visual, out = out,
          res_name = "chloroplast", ex_name = "microbeam", 
          graph = F, unit1 = "micrometer", unit2 = "min",
          consecutive = grid_par[i,1], period = grid_par[i,2], fold = grid_par[i,3])
  mvtime <- as.data.frame(data.table::fread(paste0(out, "/nomodel_chloroplast_mvtime.csv")))
  no_data_full[i] <- nrow(mvtime)
  mvtime <- mvtime[!is.na(mvtime$start_time),]
  no_data[i] <- nrow(mvtime)
  vec_cor[i] <- cor(mvtime$visual_start_time, mvtime$start_time)
  vec_mae[i] <- mean(abs(mvtime$visual_start_time - mvtime$start_time))
  vec_rmse[i] <- sqrt(mean((mvtime$visual_start_time - mvtime$start_time)^2))
  no_perfect[i] <- sum(mvtime$visual_start_time == mvtime$start_time)
  print(paste(i, "in", nrow(grid_par), sep = " "))
}

# Save output
df_cor <- data.frame(consecutive = grid_par[,1],
                     period = grid_par[,2], 
                     fold = grid_par[,3], 
                     cor = vec_cor,
                     mae = vec_mae,
                     rmse = vec_rmse,
                     no_perfect = no_perfect,
                     no_data = no_data,
                     no_data_full = no_data_full)
data.table::fwrite(df_cor, file = paste0(out, "/grid_cor_mae_rmse.csv"))

# ggplot(data = df_cor, aes(x = period, y = fold, fill = cor)) +
#   geom_tile() +
#   viridis::scale_fill_viridis(option="plasma")



### Predict movement (best parameters)
df_cor <- as.data.frame(data.table::fread(file = paste0(out, "/grid_cor_mae_rmse.csv")))
df_cor <- df_cor[df_cor$no_data>100,]
#df_cor <- df_cor[df_cor$consecutive==1,]
df_cor_best1 <- df_cor[which.max(df_cor$cor),]
df_cor_best2 <- df_cor[which.min(df_cor$mae),]
df_cor_best3 <- df_cor[which.min(df_cor$rmse),]
df_cor_best4 <- df_cor[which.max(df_cor$no_perfect),]
df_cor_best <- df_cor_best3
#list_cor_best <- list(df_cor_best1, df_cor_best2, df_cor_best3)
#df_cor_best <- list_cor_best[[which.max(c(list_cor_best[[1]]$no_perfect, 
#                                          list_cor_best[[2]]$no_perfect, 
#                                          list_cor_best[[3]]$no_perfect))]]
nomodel(cell_list = cell_list, visual = visual, out = out,
        res_name = "chloroplast", ex_name = "microbeam", 
        unit1 = "micrometer", unit2 = "min",
        consecutive = df_cor_best$consecutive, 
        period = df_cor_best$period, fold = df_cor_best$fold)
mvtime <- as.data.frame(data.table::fread(paste0(out, "/nomodel_chloroplast_mvtime.csv")))
mvtime <- mvtime[!is.na(mvtime$start_time),]
cor(mvtime$visual_start_time, mvtime$start_time)
plot(mvtime$visual_start_time, mvtime$start_time)
abline(0,1)
nrow(mvtime)

mvtime <- mvtime[1:max(which(mvtime$cell==4)),]
cor(mvtime$visual_start_time, mvtime$start_time)
plot(mvtime$visual_start_time, mvtime$start_time)

