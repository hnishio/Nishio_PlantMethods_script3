
# Load packages
library(cellssm)
library(tictoc)
cmdstanr::set_cmdstan_path("~/cmdstan/")

# Create output directory
out <- "06_ssm_individual_calctime"
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

# Execution of state-space modeling
tic()
ssm_individual(cell_list = cell_list, out = out,
               res_name = "chloroplast", ex_name = "microbeam", 
               unit1 = "micrometer", unit2 = "min",
               stepwise = c(99, 90), graph = F, diagnosis = F)
tictoc_time <- toc()
calc_time <- data.frame(
  calc_time = as.numeric(tictoc_time$toc - tictoc_time$tic))
data.table::fwrite(calc_time, file = paste0(out, "/ssm_individual_calctime_cell2.csv"))

