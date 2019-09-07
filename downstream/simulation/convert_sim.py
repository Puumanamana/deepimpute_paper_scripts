import pandas as pd
import numpy as np

raw = pd.read_csv("raw.csv",index_col=0).T

np.save("raw.npy",raw)
np.save("cells.npy",raw.index)
np.save("genes.npy",raw.columns)

metadata = pd.read_csv("cell_info.csv",index_col=0)
metadata = pd.DataFrame(metadata.Group)
metadata.columns = ["celltype"]
metadata.to_csv("metadata_noNA.csv")

print("Data ready")
