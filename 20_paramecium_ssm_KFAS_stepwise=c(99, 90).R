
# Load packages
library(cellssm)

# Create output directory
out <- "20_paramecium_ssm_KFAS_stepwise=c(99, 90)"
if(file.exists(out)==F){
  dir.create(out, recursive=T)
}

# Load data
data("Paramecium")
cell_list <- list(Paramecium)

# Execution of state-space modeling
ssm_KFAS(cell_list = cell_list, out = out,
         ex_sign = "positive", df_name = "experiment",
         res_name = "Paramecium", ex_name = "heat",
         unit1 = "millimeter", unit2 = "sec",
         stepwise = c(99, 90))

