library(splatter)

# params <- splatEstimate(read.csv("../Hrvatin/raw.csv",row.names=1))

# Choose a set of custom parameters
params.custom <- list(
    nGenes=4000, batchCells=2000,
    group.prob=c(.1,.1,.2,.2,.4),
    dropout.shape=-0.5, dropout.mid=1, # These parameters are described in the supplementary figures of the Splatter package. It corresponds to the shape and mean value of the logistic dropout function
    dropout.type="experiment"    
)

params <- newSplatParams()
params <- setParams(params, update=params.custom)

sim <- splatSimulateGroups(params)
matrices <- assays(sim)

raw <- matrices[["counts"]]
truth <- matrices[["TrueCounts"]]
dropout <- matrices[["Dropout"]]

# Write metadata
write.csv(colData(sim),"cell_info.csv")
write.csv(rowData(sim),"gene_info.csv")

# Write data
write.csv(raw,"raw.csv")
write.csv(truth,"imputed/truth_raw.csv")
write.csv(dropout,"dropout_mask.csv")

# Save global R object just in case
saveRDS(sim,"splatterSim.RDS")

system("python convert_sim.py")

system("/home/carisdak/downstream/scanpy/run_all.sh sim")
