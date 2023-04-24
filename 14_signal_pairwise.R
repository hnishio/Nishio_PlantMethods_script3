
# Load packages
library(cellssm)
library(tidyverse)
library(patchwork)
library(gt)

# Create output directory
out <- "14_signal_pairwise"
if(file.exists(out)==F){
  dir.create(out, recursive=T)
}

# Load plot function
source("functions/Plot_functions.R")

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
mvtime <- as.data.frame(data.table::fread("04_nomodel/nomodel_chloroplast_mvtime.csv"))


## Edit data
mvtime$distance <- as.numeric(c(cell1[30,-(1:2)], cell2[30,-(1:2)], cell3[30,-(1:2)], cell4[30,-(1:2)],
                                cell5[30,-(1:2)], cell6[30,-(1:2)], cell7[30,-(1:2)], cell8[30,-(1:2)],
                                cell9[30,-(1:2)]))
mvtime <- mvtime[!is.na(rowSums(mvtime)),]


## Line plot pairwise (x: start time, y: distance)
g_pairwise1 <- pairwise_start_dis(df=mvtime, i=1)
g_pairwise2 <- pairwise_start_dis(df=mvtime, i=2)
g_pairwise3 <- pairwise_start_dis(df=mvtime, i=3)
g_pairwise4 <- pairwise_start_dis(df=mvtime, i=4)
g_pairwise5 <- pairwise_start_dis(df=mvtime, i=5)
g_pairwise6 <- pairwise_start_dis(df=mvtime, i=6)
g_pairwise7 <- pairwise_start_dis(df=mvtime, i=7)
g_pairwise8 <- pairwise_start_dis(df=mvtime, i=8)
g_pairwise9 <- pairwise_start_dis(df=mvtime, i=9)


## Line plot (x: start time, y: distance) for aligned chloroplasts
df_align <- as.data.frame(data.table::fread("data/chl_align.csv"))

g_align1 <- aligned_start_dis(df=mvtime, df_align=df_align, i=1)
g_align2 <- aligned_start_dis(df=mvtime, df_align=df_align, i=2)
g_align3 <- aligned_start_dis(df=mvtime, df_align=df_align, i=3)
g_align4 <- aligned_start_dis(df=mvtime, df_align=df_align, i=4)
g_align5 <- aligned_start_dis(df=mvtime, df_align=df_align, i=5)
g_align6 <- aligned_start_dis(df=mvtime, df_align=df_align, i=6)
g_align7 <- aligned_start_dis(df=mvtime, df_align=df_align, i=7)
g_align8 <- aligned_start_dis(df=mvtime, df_align=df_align, i=8)
g_align9 <- aligned_start_dis(df=mvtime, df_align=df_align, i=9)


## Integration of all plots into a figure
void <- ggplot() + theme_void() + theme(plot.tag = element_text(size = 12, face = "bold"))

g <- (void + labs(tag = "A")) - 
  ((g_align1 + labs(tag = "B")) + g_align2 + g_align3 +
     g_align4 + g_align5 + g_align6 + g_align7 + g_align8 + g_align9 + 
     (g_pairwise1 + labs(tag = "C")) + g_pairwise2 + g_pairwise3 + g_pairwise4 +
     g_pairwise5 + g_pairwise6 + g_pairwise7 + g_pairwise8 + g_pairwise9 +
     plot_layout(ncol = 3)) + 
  plot_layout(widths = c(0.5, 1))

ggsave(paste0(out, "/signal_pairwise.pdf"),
       g, width = 160, height = 247, units = "mm")



### Chi square test

## All combinations
no_combs <- NULL
no_posis <- NULL
no_negas <- NULL

for(i in 1:9){
  data <- mvtime[mvtime$cell==i,]
  align <- as.data.frame(t(combn(x = data$chloroplast, m = 2)))
  names(align) <- c("from", "to")
  no_combs <- c(no_combs, nrow(align))
  
  no_positive <- 0
  no_negative <- 0
  for(j in 1:nrow(align)){
    x_from <- data[data$chloroplast == align$from[j],"start_time"]
    x_to <- data[data$chloroplast == align$to[j],"start_time"]
    y_from <- data[data$chloroplast == align$from[j],"distance"]
    y_to <- data[data$chloroplast == align$to[j],"distance"]
    if((y_to - y_from) / (x_to - x_from) >= 0){
      no_positive <- no_positive + 1
    }else{
      no_negative <- no_negative + 1
    }
  }#for(j in 1:nrow(align)){
  no_posis <- c(no_posis, no_positive)
  no_negas <- c(no_negas, no_negative)
  
}#for(i in 1:9){

sum_no_combs <- sum(no_combs)
sum_no_posis <- sum(no_posis)
sum_no_negas <- sum(no_negas)



## aligned
no_combs <- NULL
no_posis <- NULL
no_negas <- NULL

for(i in 1:9){
  data <- mvtime[mvtime$cell==i,]
  align <- df_align[df_align$cell==i,]
  
  rem_chl <- base::setdiff(unique(c(align$from, align$to)), data$chl)
  if(length(rem_chl) != 0){
    align <- align[!(align$from==rem_chl | align$to==rem_chl),]
  }
  no_combs <- c(no_combs, nrow(align))
  
  no_positive <- 0
  no_negative <- 0
  for(j in 1:nrow(align)){
    if(sum(data$chloroplast == align$from[j]) > 0 & sum(data$chloroplast == align$to[j]) > 0){
      x_from <- data[data$chloroplast == align$from[j],"start_time"]
      x_to <- data[data$chloroplast == align$to[j],"start_time"]
      y_from <- data[data$chloroplast == align$from[j],"distance"]
      y_to <- data[data$chloroplast == align$to[j],"distance"]
      if((y_to - y_from) / (x_to - x_from) >= 0){
        no_positive <- no_positive + 1
      }else{
        no_negative <- no_negative + 1
      }
    }#if(sum(data$chloroplast ==
  }#for(j in 1:nrow(align)){
  no_posis <- c(no_posis, no_positive)
  no_negas <- c(no_negas, no_negative)
  
}#for(i in 1:9){

sum_no_align <- sum(no_combs)
sum_no_align_posis <- sum(no_posis)
sum_no_align_negas <- sum(no_negas)

sum_no_notalign_posis <- sum_no_posis - sum_no_align_posis
sum_no_notalign_negas <- sum_no_negas - sum_no_align_negas
sum_no_notalign <- sum_no_notalign_posis + sum_no_notalign_negas

#                               sign
#                               posi                    nega  
# align   yes     sum_no_align_posis      sum_no_align_negas   sum_no_align
#         no   sum_no_notalign_posis   sum_no_notalign_negas   sum_no_notalign
#                       sum_no_posis            sum_no_negas   sum_no_combs

df_chi <- data.frame(align = rep(c("yes", "no"), 2),
                     sign = rep(c("posi", "nega"), each = 2),
                     num = c(sum_no_align_posis, sum_no_notalign_posis,
                             sum_no_align_negas, sum_no_notalign_negas))
xtab_chi <- stats::xtabs(formula = num~., data = df_chi)

(fisher_res <- fisher.test(xtab_chi))
names(fisher_res)



### Visualization of cross-table
df_tab <- data.frame(align = c("Yes", "No"),
                     posi = c(sum_no_align_posis, sum_no_notalign_posis),
                     nega = c(sum_no_align_negas, sum_no_notalign_negas))

gt_tab <- df_tab %>% 
  gt() %>%
  cols_label(align = md("**Align**"), posi = md("Positive"), nega = md("Negative"))  %>% # column names
  tab_spanner(label = md("**Slope**"), columns = c("posi", "nega")) %>%
  tab_spanner(label = "", columns = c("align")) %>%
  tab_options(table.width = pct(0),
              data_row.padding = px(4),
              table_body.hlines.width = 0,
              column_labels.border.top.width = 0,
              column_labels.border.bottom.width = 2,
              table.border.top.width = 0,
              table.border.bottom.width = 0,
              column_labels.border.top.color = "black",
              column_labels.border.bottom.color = "black",
              table.border.top.color = "black",
              table_body.border.top.color = "black",
              table_body.border.bottom.color = "black",
              table.font.size = px(9.33),
              column_labels.font.size = px(9.33),
              source_notes.font.size = px(9.33)) %>% 
  opt_table_font(
    font = c("Helvetica")) %>%
  cols_align(align = "center", 
             columns = c("posi", "nega")) %>%
  tab_source_note(source_note = paste(md("P"), " = ", round(fisher_res$p.value, 4), 
                                      " (Fisher's exact test)", sep=""))

gtsave(gt_tab, paste0(out, "/crosstable.pdf"))

