import pandas as pd
import numpy as np
import h5py

import argparse
from scipy.stats import pearsonr

from deepimpute.multinet import MultiNet

import os,sys
PARENT_DIR = os.path.join(sys.path[0], '..')
sys.path.insert(1, PARENT_DIR)

#------------------------# Parse args #------------------------#

parser = argparse.ArgumentParser(description='Impute data.')
parser.add_argument('-d', type=str, default='neuron9k')

args = parser.parse_args()

dataset = args.d
ncores = 40
n_iter = 10

#------------------------# Data #------------------------#

handle = h5py.File('{}/paper_data/accuracy.h5'.format(PARENT_DIR),'r').get(dataset)

cells = handle.get('cells')[:].astype(str)
genes = handle.get('genes')[:].astype(str)

raw = pd.DataFrame(handle.get('raw')[:], index=cells, columns=genes)
truth = pd.DataFrame(handle.get('truth')[:], index=cells, columns=genes)

mask = (raw != truth)

print('Raw data loaded')

outputdir = "{}/results/training_w_subsets".format(PARENT_DIR)
if not os.path.exists(outputdir):
    os.mkdir(outputdir)

#-------------------------# DeepImpute #-------------------------#

cellRatios = [0.05,0.1,0.2,0.4,0.6,0.8,0.9,1]

for cellRatio in cellRatios:
    for nb in range(n_iter):
        print("Cellratio: {} (iteration {})".format(cellRatio,nb))

        output_name = "{}/deepImpute_{}_{}.npy".format(outputdir,cellRatio,nb)
        
        if not os.path.exists(output_name):
            model = MultiNet(ncores=ncores,verbose=0)
            model.fit(raw,cell_subset=cellRatio)
            imputed = model.predict(raw,imputed_only=True)
            np.save(output_name,imputed.values)
            np.save('{}/imputed_genes.npy'.format(outputdir),imputed.columns)

imputed_genes = np.load('{}/imputed_genes.npy'.format(outputdir))
truth = truth[imputed_genes].values
mask = mask[imputed_genes].values

#-------------------------# MSE and Pearson #-------------------------#

scores = []

for ratio in cellRatios:
    for nb in range(n_iter):
        imputed = np.load("{}/deepImpute_{}_{}.npy".format(outputdir,ratio,nb))

        x_true = np.log1p(truth[mask])
        x_imputed = np.log1p(imputed[mask])

        pears = pearsonr(x_true,x_imputed)[0]
        MSE = np.mean((x_true-x_imputed)**2)

        info = [ ratio, nb, pears, MSE ]
        scores.append(info)

        print('{0:>10}: Pearson = {1:1.3f}, MSE = {2:.3f}'.
              format(ratio,pears,MSE))

scores = pd.DataFrame(scores, columns=['ratio','iter','Pearson','MSE'])
scores = scores.melt(id_vars=['ratio','iter'])
scores.columns = ['ratio','iter','metric','score']

scores.to_csv('{}/metrics_with_increasing_fraction.csv'.format(outputdir))
