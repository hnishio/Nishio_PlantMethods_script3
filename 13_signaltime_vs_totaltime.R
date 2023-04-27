
# Load packages
library(cellssm)
library(ggplot2)
library(patchwork)

# Create output directory
out <- "13_signaltime_vs_totaltime"
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


## Edit data
mvtime$distance <- as.numeric(c(cell1[30,-(1:2)], cell2[30,-(1:2)], cell3[30,-(1:2)], cell4[30,-(1:2)],
                                cell5[30,-(1:2)], cell6[30,-(1:2)], cell7[30,-(1:2)], cell8[30,-(1:2)],
                                cell9[30,-(1:2)]))
mvtime <- mvtime[!is.na(rowSums(mvtime)),]


## Linear regression
for(i in 1:9){
  eval(parse(text = paste0("
    lm_cell", i, " <- RobustLinearReg::siegel_regression(distance ~ start_time, data = mvtime[mvtime$cell==", i, ",])"
  )))
}
lm_allcells <- RobustLinearReg::siegel_regression(distance ~ start_time, data = mvtime)

medsig <- median(c(lm_cell1$coefficients[2], lm_cell2$coefficients[2], lm_cell3$coefficients[2], 
                   lm_cell4$coefficients[2], lm_cell5$coefficients[2], lm_cell6$coefficients[2], 
                   lm_cell7$coefficients[2], lm_cell8$coefficients[2], lm_cell9$coefficients[2]))
## medsig: 0.7716529

## Signal transfer time vs. Total reaction time
# for(i in 1:9){
#   eval(parse(text = paste0("
#     signaling_time_cell", i, " <- mvtime[mvtime$cell==", i, ",]$distance / lm_cell", i, "$coefficients[2]
#     start_time_cell", i, " <- mvtime[mvtime$cell==", i, ",]$start_time"
#   )))
# }
# signaling_time_allcells <- mvtime$distance / lm_allcells$coefficients[2]
# start_time_allcells <- mvtime$start_time

for(i in 1:9){
  eval(parse(text = paste0("
    signaling_time_cell", i, " <- mvtime[mvtime$cell==", i, ",]$distance / medsig
    start_time_cell", i, " <- mvtime[mvtime$cell==", i, ",]$start_time"
  )))
}
signaling_time_allcells <- mvtime$distance / medsig
start_time_allcells <- mvtime$start_time


## Boxplot
df_cat_signaling <- data.frame(
  cell=c(rep("Cell 1", length(start_time_cell1)),
         rep("Cell 2", length(start_time_cell2)),
         rep("Cell 3", length(start_time_cell3)), 
         rep("Cell 4", length(start_time_cell4)),
         rep("Cell 5", length(start_time_cell5)),
         rep("Cell 6", length(start_time_cell6)),
         rep("Cell 7", length(start_time_cell7)), 
         rep("Cell 8", length(start_time_cell8)),
         rep("Cell 9", length(start_time_cell9)),
         rep("Cell 1", length(signaling_time_cell1)),
         rep("Cell 2", length(signaling_time_cell2)),
         rep("Cell 3", length(signaling_time_cell3)), 
         rep("Cell 4", length(signaling_time_cell4)),
         rep("Cell 5", length(signaling_time_cell5)),
         rep("Cell 6", length(signaling_time_cell6)),
         rep("Cell 7", length(signaling_time_cell7)), 
         rep("Cell 8", length(signaling_time_cell8)),
         rep("Cell 9", length(signaling_time_cell9))),
  STorSIG = c(rep("Total reaction time", 
                  length(start_time_cell1)+length(start_time_cell2)+length(start_time_cell3)+
                    length(start_time_cell4)+length(start_time_cell5)+length(start_time_cell6)+
                    length(start_time_cell7)+length(start_time_cell8)+length(start_time_cell9)), 
              rep("Signal transfer time", 
                  length(signaling_time_cell1)+length(signaling_time_cell2)+length(signaling_time_cell3)+
                    length(signaling_time_cell4)+length(signaling_time_cell5)+length(signaling_time_cell6)+
                    length(signaling_time_cell7)+length(signaling_time_cell8)+length(signaling_time_cell9))),
  time=c(start_time_cell1, start_time_cell2, start_time_cell3, 
         start_time_cell4, start_time_cell5, start_time_cell6,
         start_time_cell7, start_time_cell8, start_time_cell9, 
         signaling_time_cell1, signaling_time_cell2, signaling_time_cell3, 
         signaling_time_cell4, signaling_time_cell5, signaling_time_cell6,
         signaling_time_cell7, signaling_time_cell8, signaling_time_cell9))
df_cat_signaling$STorSIG <- factor(df_cat_signaling$STorSIG, levels=c("Signal transfer time", "Total reaction time"))

## Wilcoxon test
#wilcox.exact(signaling_time_cell2, start_time_cell2, "less", paired = T)
#wilcox.exact(signaling_time_cell3, start_time_cell3, "less", paired = T)
#wilcox.exact(signaling_time_cell4, start_time_cell4, "less", paired = T)

g_sig_st <- ggplot(df_cat_signaling, aes(x = cell, y = time, fill = STorSIG)) +
  geom_boxplot(outlier.shape = NA) +
  #geom_jitter(size = 0.02, alpha = 0.3) +
  geom_point(size = 0.1, alpha = 0.5, position=position_jitterdodge()) +
  scale_fill_manual(values=c("darkturquoise", "gray70"))+
  coord_cartesian(ylim = c(0, 57), clip = "off") +
  annotate("text", x = 5, y = 57, label="Hypothesis: signal transfer time < total reaction time", size = 7/ggplot2::.pt) +
  annotate("text", x = 1:9, y = rep(51.5, 9), label=rep("NS",9), size = 7/ggplot2::.pt) +
  annotate("segment", x = 0.7, xend = 1.3, y = 49, yend = 49, size = 0.5,
           arrow = arrow(ends = "both", length = unit(0, "mm"))) +
  annotate("segment", x = 1.7, xend = 2.3, y = 49, yend = 49, size = 0.5,
           arrow = arrow(ends = "both", length = unit(0, "mm"))) +
  annotate("segment", x = 2.7, xend = 3.3, y = 49, yend = 49, size = 0.5,
           arrow = arrow(ends = "both", length = unit(0, "mm"))) +
  annotate("segment", x = 3.7, xend = 4.3, y = 49, yend = 49, size = 0.5,
           arrow = arrow(ends = "both", length = unit(0, "mm"))) +
  annotate("segment", x = 4.7, xend = 5.3, y = 49, yend = 49, size = 0.5,
           arrow = arrow(ends = "both", length = unit(0, "mm"))) +
  annotate("segment", x = 5.7, xend = 6.3, y = 49, yend = 49, size = 0.5,
           arrow = arrow(ends = "both", length = unit(0, "mm"))) +
  annotate("segment", x = 6.7, xend = 7.3, y = 49, yend = 49, size = 0.5,
           arrow = arrow(ends = "both", length = unit(0, "mm"))) +
  annotate("segment", x = 7.7, xend = 8.3, y = 49, yend = 49, size = 0.5,
           arrow = arrow(ends = "both", length = unit(0, "mm"))) +
  annotate("segment", x = 8.7, xend = 9.3, y = 49, yend = 49, size = 0.5,
           arrow = arrow(ends = "both", length = unit(0, "mm"))) +
  theme_bw() +
  theme(plot.title=element_text(size=7),
        axis.title.y=element_text(size=7), 
        axis.title.x=element_blank(),
        axis.text=element_text(size=7),
        #axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = c(0.77, 0.7),
        legend.background = element_rect(fill="transparent"),
        legend.title = element_blank(),
        legend.text = element_text(size=7),
        #legend.margin=margin(0, 0, 0, -5),
        legend.spacing.x = unit(0.7, 'mm'),
        legend.key.width = unit(2, "mm"),
        plot.tag = element_text(size = 12, face = "bold")) + 
  labs(y="Time (min)")



##### Explanation of reaction time
g_text <- ggplot() + theme_void() + 
  geom_segment(aes(x = -1, y = 0, xend = 1.05, yend = 0),
               arrow = arrow(length = unit(0.3, "cm")), size = 1) +
  geom_curve(aes(x = -1, y = 0, xend = -0.7, yend = 0), curvature = -0.6, size = 0.2) +
  geom_curve(aes(x = -0.7, y = 0, xend = 0.7, yend = 0), curvature = -0.2, size = 0.2) +
  geom_curve(aes(x = 0.7, y = 0, xend = 1, yend = 0), curvature = -0.6, size = 0.2) +
  geom_curve(aes(x = -1, y = 0, xend = 1, yend = 0), curvature = 0.2, size = 0.2) +
  geom_segment(aes(x = -1, y = 0.2, xend = -1, yend = -0.6), lty = 2, size = 0.2) +
  geom_segment(aes(x = 1.05, y = 0.2, xend = 1.05, yend = -0.6), lty = 2, size = 0.2) +
  theme(plot.tag = element_text(size = 12, face = "bold")) +
  coord_cartesian(xlim = c(-1, 1.3), ylim = c(-1.5, 0.5), clip = "off") +
  annotate("text", x=0, y=-0.44, label = "Total reaction time\n(= start time)", size = 7/ggplot2::.pt) +
  annotate("text", x=-0.85, y=0.3, label = "Warm-up\ntime", size = 7/ggplot2::.pt) + 
  annotate("text", x=0, y=0.3, label = "Signal transfer time", size = 7/ggplot2::.pt) + 
  annotate("text", x=0.85, y=0.3, label = "Warm-up\ntime", size = 7/ggplot2::.pt) +
  annotate("text", x=-1, y=-0.72, label = "Start of \nmicrobeam", size = 7/ggplot2::.pt) + 
  annotate("text", x=1.05, y=-0.79, label = "Start of \nchloroplast\naccumulation", size = 7/ggplot2::.pt)


## Integrate figures
g <- (g_text + labs(tag = "A")) + (g_sig_st + labs(tag = "B")) +
  plot_layout(ncol = 2, widths = c(1,1.3))
set.seed(1)
suppressWarnings(ggsave(paste0(out, "/signaltime_vs_totaltime.pdf"),
                        g, height = 70, width = 183, units = "mm"))

