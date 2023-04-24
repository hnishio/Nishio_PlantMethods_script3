
# Load packages
library(cellssm)
library(ggplot2)
library(patchwork)

# Create output directory
out <- "11_lm_signal_nomodel"
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





##### robust = T #####

## Edit data
mvtime <- as.data.frame(data.table::fread("04_nomodel/nomodel_chloroplast_mvtime.csv"))
mvtime$distance <- as.numeric(c(cell1[30,-(1:2)], cell2[30,-(1:2)], cell3[30,-(1:2)], cell4[30,-(1:2)],
                                cell5[30,-(1:2)], cell6[30,-(1:2)], cell7[30,-(1:2)], cell8[30,-(1:2)],
                                cell9[30,-(1:2)]))
mvtime <- mvtime[!is.na(rowSums(mvtime)),]

## Boxplot of signaling speed
for(i in 1:9){
  eval(parse(text = paste0("
    lm_cell", i, " <- RobustLinearReg::siegel_regression(distance ~ start_time, data = mvtime[mvtime$cell==", i, ",])"
  )))
}
signalingspeed <- c(lm_cell1$coefficients[2], lm_cell2$coefficients[2], lm_cell3$coefficients[2], 
                    lm_cell4$coefficients[2], lm_cell5$coefficients[2], lm_cell6$coefficients[2], 
                    lm_cell7$coefficients[2], lm_cell8$coefficients[2], lm_cell9$coefficients[2])
df_sig <- data.frame(category = rep("Cell 1-9", 9),
                     signalingspeed = signalingspeed)

g_sig <- ggplot(df_sig, aes(x = category, y = signalingspeed)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(size = 0.8, alpha = 0.5) +
  theme_bw() +
  theme(plot.title=element_text(size=7),
        axis.title.y=element_text(size=7), 
        axis.title.x=element_blank(),
        axis.text=element_text(size=7),
        plot.tag = element_text(size = 12, face = "bold")) + 
  labs(y=bquote(paste("Signal transfer speed ", (mu*m/min), sep = "")))


## nomodel
mvtime <- as.data.frame(data.table::fread("04_nomodel/nomodel_chloroplast_mvtime.csv"))
glist <- lm_signal(cell_list = cell_list, mvtime = mvtime, robust = T,
                   ex_name = "microbeam", unit1 = "micrometer", unit2 = "min")

void <- ggplot() + theme_void()
g <- (glist[[1]] + labs(tag = "A")) + (glist[[2]] + labs(tag = "B")) + (glist[[3]] + labs(tag = "C")) + 
  (glist[[4]] + labs(tag = "D")) + (glist[[5]] + labs(tag = "E")) + (glist[[6]] + labs(tag = "F")) + 
  (glist[[7]] + labs(tag = "G")) + (glist[[8]] + labs(tag = "H")) + (glist[[9]] + labs(tag = "I")) + 
  (glist[[10]] + labs(tag = "J")) + (glist[[11]] + labs(tag = "K")) +
  ((g_sig + labs(tag = "L")) + void + plot_layout(ncol = 2, widths = c(1,1.5))) +
  plot_layout(ncol = 3) &
  theme(plot.tag = element_text(size = 12, face = "bold"))
set.seed(6)
suppressWarnings(ggsave(paste0(out, "/lm_signal_nomodel_robustT_consecutive==1-5_rmse.pdf"),
                        g, height = 220, width = 50*3, units = "mm"))





##### robust = F #####

## Edit data
mvtime <- as.data.frame(data.table::fread("04_nomodel/nomodel_chloroplast_mvtime.csv"))
mvtime$distance <- as.numeric(c(cell1[30,-(1:2)], cell2[30,-(1:2)], cell3[30,-(1:2)], cell4[30,-(1:2)],
                                cell5[30,-(1:2)], cell6[30,-(1:2)], cell7[30,-(1:2)], cell8[30,-(1:2)],
                                cell9[30,-(1:2)]))
mvtime <- mvtime[!is.na(rowSums(mvtime)),]

## Boxplot of signaling speed
for(i in 1:9){
  eval(parse(text = paste0("
    lm_cell", i, " <- stats::lm(distance ~ start_time, data = mvtime[mvtime$cell==", i, ",])"
  )))
}
signalingspeed <- c(lm_cell1$coefficients[2], lm_cell2$coefficients[2], lm_cell3$coefficients[2], 
                    lm_cell4$coefficients[2], lm_cell5$coefficients[2], lm_cell6$coefficients[2], 
                    lm_cell7$coefficients[2], lm_cell8$coefficients[2], lm_cell9$coefficients[2])
df_sig <- data.frame(category = rep("Cell 1-9", 9),
                     signalingspeed = signalingspeed)

g_sig <- ggplot(df_sig, aes(x = category, y = signalingspeed)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(size = 0.8, alpha = 0.5) +
  theme_bw() +
  theme(plot.title=element_text(size=7),
        axis.title.y=element_text(size=7), 
        axis.title.x=element_blank(),
        axis.text=element_text(size=7),
        plot.tag = element_text(size = 12, face = "bold")) + 
  labs(y=bquote(paste("Signal transfer speed ", (mu*m/min), sep = "")))


## nomodel
mvtime <- as.data.frame(data.table::fread("04_nomodel/nomodel_chloroplast_mvtime.csv"))
glist <- lm_signal(cell_list = cell_list, mvtime = mvtime, robust = F,
                   ex_name = "microbeam", unit1 = "micrometer", unit2 = "min")

void <- ggplot() + theme_void()
g <- (glist[[1]] + labs(tag = "A")) + (glist[[2]] + labs(tag = "B")) + (glist[[3]] + labs(tag = "C")) + 
  (glist[[4]] + labs(tag = "D")) + (glist[[5]] + labs(tag = "E")) + (glist[[6]] + labs(tag = "F")) + 
  (glist[[7]] + labs(tag = "G")) + (glist[[8]] + labs(tag = "H")) + (glist[[9]] + labs(tag = "I")) + 
  (glist[[10]] + labs(tag = "J")) + (glist[[11]] + labs(tag = "K")) + 
  ((g_sig + labs(tag = "L")) + void + plot_layout(ncol = 2, widths = c(1,1.5))) +
  plot_layout(ncol = 3) &
  theme(plot.tag = element_text(size = 12, face = "bold"))
set.seed(6)
suppressWarnings(ggsave(paste0(out, "/lm_signal_nomodel_robustF_consecutive==1-5_rmse.pdf"),
                        g, height = 220, width = 50*3, units = "mm"))

