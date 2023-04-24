
# Load packages
library(cellssm)

# Create output directory
out <- "17_nucleus_ssm_KFAS_stepwise=c(99, 90)"
if(file.exists(out)==F){
  dir.create(out, recursive=T)
}

# Load data
cell1 <- read.csv("data/nucleus1.csv")
cell2 <- read.csv("data/nucleus2.csv")
cell_list <- list(cell1, cell2)

# Execution of state-space modeling
ssm_KFAS(cell_list = cell_list, out = out,
         res_name = "nucleus", ex_name = "microbeam", 
         unit1 = "micrometer", unit2 = "min",
         stepwise = c(99, 90))

