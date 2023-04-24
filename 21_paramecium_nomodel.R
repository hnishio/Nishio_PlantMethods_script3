
# Load packages
library(cellssm)

# Create output directory
out <- "21_paramecium_nomodel"
if(file.exists(out)==F){
  dir.create(out, recursive=T)
}

# Load data
data("Paramecium")
cell_list <- list(Paramecium)

# Execution of state-space modeling
nomodel(cell_list = cell_list, out = out,
        ex_sign = "positive", fold = 1, df_name = "experiment",
        res_name = "Paramecium", ex_name = "heat",
        unit1 = "millimeter", unit2 = "sec")

