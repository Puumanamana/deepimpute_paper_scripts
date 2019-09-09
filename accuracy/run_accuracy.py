import os,sys
import pandas as pd
import numpy as np
import argparse
import h5py
import scipy.stats

from deepimpute.multinet import MultiNet

PARENT_DIR = os.path.join(sys.path[0], '..')
sys.path.insert(1, PARENT_DIR)

from config import imputation_methods

#------------------------# Parse args #------------------------#

parser = argparse.ArgumentParser(description='Impute data.')
parser.add_argument('-d', type=str, default='jurkat')
args = parser.parse_args()

dataset = args.d

#------------------------# raw and true data #------------------------#

handle = h5py.File('{}/paper_data/accuracy.h5'.format(PARENT_DIR),'r').get(dataset)

cells = handle.get('cells')[:].astype(str)
genes = handle.get('genes')[:].astype(str)

raw = pd.DataFrame(handle.get('raw'), index=cells, columns=genes)
truth = pd.DataFrame(handle.get('truth'), index=cells, columns=genes)

mask = (raw != truth)

#-------------------------# DeepImpute #-------------------------#

model = MultiNet(ncores=20)
model.fit(raw)

df_deepImpute = model.predict(raw,imputed_only=True)
gene_subset = df_deepImpute.columns

imputation = { 'deepImpute': np.log1p(df_deepImpute) }

#------------------------# Import other imputation matrices #------------------------#

for method in imputation_methods[1:]:
    print("Loading {}".format(method))
    imputed = pd.DataFrame(handle.get(method), index=cells, columns=genes)
    imputation[method] = np.log1p(imputed[gene_subset])

truth = np.log1p(truth.reindex(columns=gene_subset))
raw = np.log1p(raw.reindex(columns=gene_subset))
mask = mask.reindex(columns=gene_subset)

METHODS = list(imputation.keys())
print("Loaded: ", METHODS)

#------------------------# Global metrics #------------------------#

if not os.path.exists("{}/results/accuracy".format(PARENT_DIR)):
    os.mkdir("{}/results/accuracy".format(PARENT_DIR))

scatter_data = { m: np.log1p(imputation[m].values[mask.values])
                 for m in METHODS }
scatter_data["truth"] = np.log1p(truth.values[mask.values])

scatter_data = pd.melt(pd.DataFrame(scatter_data), id_vars=['truth'])
scatter_data['dataset'] = [dataset] * scatter_data.shape[0]
scatter_data.columns = ["Truth","method","Imputed","dataset"]
scatter_data.to_csv("{}/results/accuracy/scatter_{}.csv".format(PARENT_DIR,dataset))

#------------------------# Metric per cell #------------------------#

MSE_cells = pd.DataFrame({ method: ((imputation[method][mask]
                                     -truth[mask])**2).mean(axis=1)
                           for method in METHODS }
                         ).dropna()

MSE_genes = pd.DataFrame({ method: ((imputation[method][mask]
                                     -truth[mask])**2).mean(axis=0)
                           for method in METHODS }
                         ).dropna()

## Significance with paired t-test
cell_pvals = pd.DataFrame([ (dataset, method, "cell",
                scipy.stats.ttest_rel(MSE_cells['deepImpute'],
                                      MSE_cells[method])[1])
                            for method in imputation_methods[1:] ],
                          columns=["dataset","method","axis","value"])
                          
gene_pvals = pd.DataFrame([ (dataset, method, "gene",
                scipy.stats.ttest_rel(MSE_genes['deepImpute'],
                                      MSE_genes[method])[1])
                            for method in imputation_methods[1:] ],
                          columns=["dataset","method","axis","value"])

output_file = "{}/results/accuracy/pvalues_table.csv".format(PARENT_DIR)
if os.path.exists(output_file):
    statistics = pd.read_csv(output_file,index_col=0)
    statistics = pd.concat([statistics,cell_pvals])
    statistics = pd.concat([statistics,gene_pvals])

else:
    statistics = pd.concat([cell_pvals,gene_pvals])

statistics.to_csv(output_file)

#------------------------# Metric per gene #------------------------#

MSE_cells["axis"] = ["cells"] * MSE_cells.shape[0]
MSE_genes["axis"] = ["genes"] * MSE_genes.shape[0]

MSE = pd.concat([pd.melt(MSE_cells, id_vars=["axis"]),
                 pd.melt(MSE_genes, id_vars=["axis"])])
MSE["dataset"] = [dataset] * MSE.shape[0]
MSE.columns = ["axis","method","MSE","dataset"]

MSE.to_csv("{}/results/accuracy/MSE_{}.csv".format(PARENT_DIR,dataset))
