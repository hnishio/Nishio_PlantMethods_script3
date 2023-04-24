
# Load packages
library(ggplot2)
library(patchwork)
library(data.table)

# Load plot functions
source("functions/Plot_functions.R")

# Create an output directory
out <- "05_visual_vs_predicted"
if(file.exists(out)==F){
  dir.create(out, recursive=T)
}


## Bayesian inference
df <- as.data.frame(data.table::fread("02_ssm_individual_stepwise=c(99, 90)/csv/ssm_individual_chloroplast_mvtime.csv"))
df <- df[,1:4]
names(df)[1:4] <- c("cell", "index", "visual", "predicted")
df <- df[!is.infinite(rowSums(df)),]
g_Bayes <- cor_vis_pred(df, "Bayesian inference")


## Kalman filter
df <- as.data.frame(data.table::fread("03_ssm_KFAS_stepwise=c(99, 90)/csv/ssm_KFAS_chloroplast_mvtime.csv"))
df <- df[,1:4]
names(df)[1:4] <- c("cell", "index", "visual", "predicted")
df <- df[!is.infinite(rowSums(df)),]
g_KF <- cor_vis_pred(df, "Kalman filter")


## No model
df <- as.data.frame(data.table::fread("04_nomodel/nomodel_chloroplast_mvtime.csv"))
df <- df[,1:4]
names(df)[1:4] <- c("cell", "index", "visual", "predicted")
df <- df[!is.na(rowSums(df)),]
g_NM <- cor_vis_pred(df, "Without model")

g <- g_Bayes + g_KF + g_NM +
  plot_layout(ncol = 3)
ggsave(paste0(out, "/visual_vs_predicted_cell1-9_consecutive==1-5_rmse.pdf"),
       g, height = 50, width = 125, units = "mm")

