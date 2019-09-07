import os,sys

import pandas as pd
import numpy as np
import sklearn.metrics

import scanpy.api as sc
import scanpy.api.pp as pp
import scanpy.api.tl as tl

from deepimpute.multinet import MultiNet

def seurat(adata,filename,dataset=None):

    output_path = "../results/downstream/UMAP/{}".format(dataset)

    if not os.path.exists(output_path):
        os.system("mkdir -p {}".format(output_path))

    pp.recipe_seurat(adata)
    
    adata.X[np.isnan(adata.X)] = 0
    adata.X[np.isinf(adata.X)] = 0    

    tl.pca(adata)
    pp.neighbors(adata, use_rep='X_pca')

    tl.umap(adata)

    np.save("{}/{}.npy".format(output_path,filename),adata.obsm['X_umap'])

    print("Seurat finished")
    return adata

def cluster(adata):
    tl.leiden(adata)
    return adata

def evaluate(adata,
             imputation_name,
             output):
    truth = adata.obs['celltype']
    pred = adata.obs['leiden']
    X = adata.obsm["X_umap"]

    scores = {'adjusted_rand_score': sklearn.metrics.adjusted_rand_score(truth,pred),
              'adjusted_mutual_info_score': sklearn.metrics.adjusted_mutual_info_score(truth,pred),
              'Fowlkes-Mallows': sklearn.metrics.fowlkes_mallows_score(truth,pred),
              'silhouette_score': sklearn.metrics.silhouette_score(X,truth.tolist())}
    print(scores)
    
    if os.path.exists(output):
        scores_df = pd.read_csv(output,index_col=0)
    else:
        scores_df = pd.DataFrame(columns=list(scores.keys()))

    scores_df.loc[imputation_name] = pd.Series(scores)
    
    scores_df.index.name = "imputation"
    scores_df.to_csv(output)

def extract_DEGs(adata,nGenes):
    tl.rank_genes_groups(adata,
                         groupby='celltype',
                         n_genes=nGenes,
                         method='wilcoxon')

    DEG_names = pd.DataFrame(adata.uns['rank_genes_groups']['names'])
    DEG_adjPvals = pd.DataFrame(adata.uns['rank_genes_groups']['pvals_adj'])
    DEG_pvals = pd.DataFrame(adata.uns['rank_genes_groups']['pvals'])    
    DEG_lFCs = pd.DataFrame(adata.uns['rank_genes_groups']['logfoldchanges'])

    res = { group: pd.DataFrame(np.array( [DEG_adjPvals[group].values,
                                           DEG_pvals[group].values,
                                           DEG_lFCs[group].values] ).T,
                                index=DEG_names[group],
                                columns=['adj_pval','pval','log_FC'])
            for group in DEG_adjPvals.columns }

    return res

    
def run_DI(raw):
    df = pd.DataFrame(raw.X)
    model = MultiNet(ncores=40)
    imputed = model.fit(df).predict(df)
        
    adata = sc.AnnData(imputed.values)
    adata.obs_names = raw.obs.index
    adata.var_names = raw.var.index
    adata.obs["celltype"] = raw.obs.celltype.values

if __name__ == '__main__':
    import pickle

    nargs = len(sys.argv)

    dataset = "sim"

    if nargs>1:
        dataset = sys.argv[1]
    
    print("Dataset: {}".format(dataset))

    DATA_DIR = "paper_data/downstream"

    print("Loading raw data.")
    
    raw = sc.read_h5ad("{}/raw_{}.h5ad".format(DATA_DIR,dataset))
    run_DI(raw)
    
    print("Starting postprocessing.")
    DEGs = {}
    scores_classif = {}

    for method in ["raw","deepImpute","DCA","MAGIC","SAVER","scImpute","DrImpute","VIPER"]:
        path = "{}/{}_{}.h5ad".format(DATA_DIR,method,dataset)
        if not os.path.exists(path):
            print('{} does not exists. Skipping'.format(path))
            continue
        
        adata = sc.read_h5ad(path)
        adata = seurat(adata,method,dataset=dataset)
        adata = cluster(adata)
        evaluate(adata,method,
                 "../results/downstream/leiden_{}_clustering_scores.csv".format(dataset))
        if dataset == 'sim':
            DEGs.update({ method: extract_DEGs(adata,nGenes=500) })
        print("{} processed.".format(method))

    if dataset == 'sim':
        with open('../results/downstream/DEGs.pickle','wb') as handle:
            pickle.dump(DEGs,handle)
                
