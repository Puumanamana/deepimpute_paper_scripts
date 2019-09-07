import os
from config import imputation_methods
import pandas as pd
import numpy as np
import argparse
import h5py

import matplotlib.pyplot as plt

from deepimpute.multinet import MultiNet

#------------------------# Parse args #------------------------#

METHODS = imputation_methods[1:5]

parser = argparse.ArgumentParser(description='Impute data.')
parser.add_argument('-d', type=str, default='jurkat')
args = parser.parse_args()

dataset = args.d

#------------------------# raw and true data #------------------------#

handle = h5py.File('../paper_data/accuracy.h5','r').get(dataset)

cells = handle.get('cells')[:].astype(str)
genes = handle.get('genes')[:].astype(str)

raw = pd.DataFrame(handle.get('raw'), index=cells, columns=genes)
truth = pd.DataFrame(handle.get('truth'), index=cells, columns=genes)

mask = (raw != truth)

outputdir = "../results/dropout_effect"
if not os.path.exists(outputdir):
    os.mkdir(outputdir)

#-------------------------# DeepImpute #-------------------------#

imputed = {}

for dp_rate in np.arange(0,1,0.1):
    print("Processing dropout rate: {:.1f}".format(dp_rate))
    
    filename = "{}/deepimpute_{:.1f}.npy".format(outputdir,dp_rate)
    geneFilename = "{}/gene_subset.npy".format(outputdir,dataset)
    
    if not os.path.exists(filename):    
        architecture = [
            {"type": "dense", "neurons": 256, "activation": "relu"},
            {"type": "dropout", "rate": dp_rate},
        ]
        model = MultiNet(architecture=architecture, ncores=40)
        model.fit(raw)
        prediction = model.predict(raw,imputed_only=True)
        
        np.save(filename,prediction.values)
        np.save(geneFilename, prediction.columns)
        
    gene_subset = np.load(geneFilename)
    imputed["{0:.1g}".format(dp_rate)] = pd.DataFrame(np.log1p(np.load(filename)),
                                                      index=cells,
                                                      columns=gene_subset)
        
#------------------------# Import other data matrices #------------------------#

truth = np.log1p(truth.reindex(columns=gene_subset))
raw = np.log1p(raw.reindex(columns=gene_subset))
mask = mask.reindex(columns=gene_subset)

#------------------------# Global metrics #------------------------#

MSE_cells = pd.DataFrame({ dp: ((imputed[dp][mask]
                                 -truth[mask])**2).mean(axis=1)
                           for dp in imputed })

MSE_genes = pd.DataFrame({ dp: ((imputed[dp][mask]
                                 -truth[mask])**2).mean(axis=0)
                           for dp in imputed })

print(MSE_cells.median())
print(MSE_genes.median())

fig,ax = plt.subplots()

ax.plot(MSE_genes.columns.values.astype(float), MSE_genes.median().tolist(),
             label="Gene level MSE", marker="x")
ax.plot(MSE_cells.columns.values.astype(float), MSE_cells.median().tolist(),
             label="Cell level MSE", marker="o")
ax.set_xlabel("Dropout rate", fontsize=12)
ax.set_ylabel("MSE score", fontsize=12)

plt.legend()
plt.title("Influence of deepImpute's dropout rate\n on imputation accuracy",
          fontsize=14, fontweight="bold")
plt.savefig("{}/dropout_effect.png".format(outputdir),dpi=300)
plt.show()
