
# Load packages
library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)
library(data.table)
library(scales)

# Create an output directory
out <- "09_compare_calctime"
if(file.exists(out)==F){
  dir.create(out, recursive=T)
}

# Load data
visual <- as.data.frame(fread("data/cell2_visual_judgetime.csv"))
bayes <- as.data.frame(fread("06_ssm_individual_calctime/ssm_individual_calctime_cell2.csv"))
kfas <- as.data.frame(fread("07_ssm_KFAS_calctime/ssm_KFAS_calctime_cell2.csv"))
nomodel <- as.data.frame(fread("08_nomodel_calctime/nomodel_calctime_cell2.csv"))

df <- data.frame(visual = sum(visual$judge_time[visual$cell==2]),
                 bayes = sum(bayes$calc_time),
                 kfas = sum(kfas$calc_time),
                 nomodel = sum(nomodel$calc_time)) %>%
  pivot_longer(cols = everything(), names_to = "method", values_to = "time")
df$method <- c("Visual", "Bayesian\ninference", "Kalman\nfilter", "Without\nmodel")

cols <- c("aquamarine3", "orange", "orange", "orange")

df$method <- factor(df$method, levels = rev(df$method))
df

g <- ggplot(df, aes(x = method, y = time, fill = method, color = method)) +
  geom_bar(width = 0.7, stat = "identity") +
  annotate("text", x=4, y=df$time[1]+200, label=comma(df$time[1]), size=7/ggplot2::.pt) +
  annotate("text", x=3, y=df$time[2]+200, label=round(df$time[2],3), size=7/ggplot2::.pt) +
  annotate("text", x=2, y=df$time[3]+200, label=round(df$time[3],3), size=7/ggplot2::.pt) +
  annotate("text", x=1, y=df$time[4]+200, label=round(df$time[4],3), size=7/ggplot2::.pt) +
  coord_flip() +
  scale_color_manual(values = rep("black",4), limits = df$method) +
  scale_fill_manual(values = cols, limits = df$method) +
  scale_y_continuous(limits = c(0, 1600), labels = comma) +
  theme_bw() +
  theme(legend.position = "none",
        axis.text = element_text(size = 7, colour = "black"),
        axis.title.y = element_blank(),
        axis.title.x = element_text(size = 7, colour = "black")) +
  labs(y = "Time required to estimate the start time of\n11 chloroplasts in cell 2 (sec)")

ggsave(paste0(out, "/compare_calctime.pdf"),
       g, height = 45, width = 67, units = "mm")

