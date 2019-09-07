import argparse

import numpy as np
import pandas as pd
import h5py
from scipy.stats import logistic

from collections import Counter

parser = argparse.ArgumentParser(description='Mask data.')
parser.add_argument('-d', type=str, default='jurkat')
args = parser.parse_args()

N_MASKED_PER_CELL = 10

dataset = args.d

def estimate_dropout_prob(data):

    def log_mean(x):
        x_nz = x[x>0]

        if len(x_nz) == 0:
            return -1e5
        else:
            return np.mean(np.log(x_nz))

    # Get mean of positive values of each gene
    log_means = np.apply_along_axis(log_mean,0,truth)
    np.mean(np.log1p(data>0), axis=0)
    # Get each gene average number of zeros
    dropout = np.mean(data==0, axis=0)
    # Fit logistic function
    params = logistic.fit((log_means,dropout))

    return params

#----------- Load data -----------#

handle = h5py.File('../paper_data/accuracy.h5','r').get(dataset)

truth = np.array(handle.get('truth'))
cells = handle.get('cells')[:].astype(str)
genes = handle.get('genes')[:].astype(str)

raw = truth.copy()

params = estimate_dropout_prob(truth)

mask = []
fails = 0

#----------- Mask each gene iteratively -----------#

for cell in range(truth.shape[0]):
    nonZeroIdx = np.nonzero(truth[cell,:])[0]
    nonZeroVals = truth[cell,nonZeroIdx]

    if len(nonZeroVals) < 50:
        fails += 1
        print("Cannot mask values for only {} cells".format(len(nonZeroVals)))
        mask.append([])
        continue
    
    probs = logistic.pdf(np.log(nonZeroVals),*params)

    mask_c = np.random.choice(
        nonZeroIdx,
        N_MASKED_PER_CELL,
        p = probs / sum(probs),
        replace=False
        )

    raw[cell,mask_c] = 0

    mask.append(mask_c)

print("Counting masked values..")

print(Counter(truth[(raw != truth)]))
print(fails)

#----------- Save masked data -----------#

handle.close()
handle = h5py.File('../paper_data/accuracy.h5','a')
handle.create_dataset('{}/raw'.format(dataset),data=raw)

raw_csv = pd.DataFrame(raw,index=cells, columns=genes).T
raw_csv.to_csv("{}_raw_to_impute.csv".format(dataset))
    
