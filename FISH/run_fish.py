#!/usr/bin/env python

import pandas as pd
import numpy as np
import os
import h5py

from deepimpute.multinet import MultiNet
import sys
sys.path.append('..')
from config import imputation_methods

#------------------------# Load FISH and raw dropseq data #------------------------#

handle = h5py.File('../paper_data/FISH.h5')
fish = pd.DataFrame(handle.get('fish/data')[:],
                    index=handle.get('fish/cells')[:].astype(str),
                    columns=handle.get('fish/genes')[:].astype(str))

cells = handle.get('dropseq/cells')[:].astype(str)
genes = handle.get('dropseq/genes')[:].astype(str)
dropseq = pd.DataFrame(handle.get('dropseq/raw')[:],
                       index=cells,columns=genes)

imputed_data = {}

#-------------------------# DeepImpute #-------------------------#

model = MultiNet(ncores=40)
model.fit(dropseq)
imputed = model.predict(dropseq,imputed_only=True)

# Only compare with imputed gene subset
genes_to_extract = np.intersect1d(imputed.columns,fish.columns)

#-------------------------# Load all imputation results #-------------------------#
print('Loading datasets')
imputed_data = { method: pd.DataFrame(handle.get('imputed/{}'.format(method))[:],
                                      index=cells, columns=genes)[genes_to_extract]
                 for method in imputation_methods[1:-1] }
imputed_data['raw'] = dropseq[genes_to_extract]
imputed_data['deepImpute'] = imputed[genes_to_extract]

print('Focusing on {} genes'.format(len(genes_to_extract)))

#-----------------# Normalization using GAPH housekeeping gene #------------------#

def normalize_dropseq_by_gapdh(df):
    """ 
    Normalize dropseq data using SAVER procedure:
    https://www.nature.com/articles/s41592-018-0033-z
    """
    lim_inf, lim_sup = df['GAPDH'].quantile(0.1), df['GAPDH'].quantile(0.9)
    
    df_normed = df.loc[(df['GAPDH']>lim_inf) & (df['GAPDH']<lim_sup)]
    
    factor = df_normed['GAPDH'] / df_normed['GAPDH'].mean()
    
    return df_normed.div(factor,axis=0)

def normalize_fish_by_gapdh(df,genes_to_extract):
    """
    Same as before but need to handle NA values in the FISH dataset
    """
    df_normed = df.copy()

    lim_inf, lim_sup = df['GAPDH'].quantile(0.1), df['GAPDH'].quantile(0.9)
    GAPDH_cond = (df['GAPDH']>lim_inf) & (df['GAPDH']<lim_sup)

    for gene in genes_to_extract:
        # Calculate the normalization factor on non-NA cells
        cells_to_keep = df.index[GAPDH_cond & ~df[gene].isnull()]
        factor = df.loc[cells_to_keep,'GAPDH'] / df.loc[cells_to_keep,'GAPDH'].mean()        
        df_normed.loc[cells_to_keep,gene] = df.loc[cells_to_keep,gene].div(factor)

    return df_normed

def gini(array,eps=1e-8):
    """Calculate the Gini coefficient of a numpy array."""
    # http://www.statsdirect.com/help/default.htm#nonparametric_methods/gini.htm
    array = np.sort(array+eps)
    index = np.arange(1,array.shape[0]+1)
    n = array.shape[0]
    
    return ((np.sum((2*index-n -1)*array)) / (n*np.sum(array)))

# Normalize by GAPDH housekeeping gene
df_normed = { met: normalize_dropseq_by_gapdh(imputed).drop('GAPDH',axis=1)
              for met,imputed in imputed_data.items() }
df_normed['FISH'] = normalize_fish_by_gapdh(fish,genes_to_extract).drop('GAPDH',axis=1)

#---------------------# GINI coefficients #----------------------#

GINI = pd.DataFrame({ met: distr.agg(lambda x: gini(x.values[~np.isnan(x.values)]))
                      for met,distr in df_normed.items() })
GINI.index.name="Gene"

# Save results
if not os.path.exists("results/fish"):
    os.mkdir("results/fish")

GINI.to_csv("results/fish/gini.csv")

#-----------------# Normalized distributions #------------------#

# Normalization
def normalize_for_distribution(imputed, fish, gene, method):
    efficiency = fish[gene].mean() / imputed[gene].mean()
    df_normed = imputed * efficiency
    
    lim_inf, lim_sup = df_normed['GAPDH'].quantile(0.1), df_normed['GAPDH'].quantile(0.9)
    df_normed = df_normed.loc[(df_normed['GAPDH']>lim_inf) & (df_normed['GAPDH']<lim_sup), :]
    df_normed_g = pd.DataFrame(df_normed[gene] / (df_normed['GAPDH'] / df_normed['GAPDH'].mean()))
    
    df_normed_g["method"] = [method] * df_normed_g.shape[0]
    df_normed_g.columns = [gene,"method"]

    return df_normed_g
    
genes_distr = []
gene_names = np.setdiff1d(genes_to_extract,["GAPDH"])

for gene in gene_names:
    res = [ normalize_for_distribution(imputed, fish, gene, method)
            for method,imputed in imputed_data.items() ]
    
    fish_res = pd.DataFrame(df_normed["FISH"][gene].dropna())
    fish_res["method"] = ["FISH"] * fish_res.shape[0]

    res += [ fish_res ]
    res = pd.concat(res, axis=0)
    res["gene"] = [gene] * res.shape[0]
    res.columns = ["value","method","gene"]
    
    genes_distr.append(res)
    
## Plot the normalized distributions
distrs = pd.concat(genes_distr).to_csv("../results/fish/normalized_distributions.csv")
    
