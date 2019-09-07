import sys
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import scanpy.api as sc

dataset = sys.argv[1]
data_dir = "../results/downstream/{}".format(dataset)

imputation = ["deepImpute","DCA","VIPER","MAGIC","SAVER","scImpute","DrImpute","raw"]

raw = sc.read_h5ad('../paper_data/downstream/raw_{}.h5ad'.format(dataset))
cells = raw.obs.index
genes = raw.var.index
metadata = raw.obs.celltype.values

components_all = []
colnames = ["UMAP_1","UMAP_2"]
for method in imputation:
    try:
        data = np.load("../results/downstream/UMAP/{}/{}.npy".format(dataset,method))
        df = pd.DataFrame(data,index=cells,columns=colnames)
        df["meta"] = metadata
        df["imputation"] = method
        components_all.append(df)
    except:
        print("Could not compute {}".format(method))
    
components_all = pd.concat(components_all)

if df.shape[0] > 5000:
    size = 1
else:
    size = 10

wrap = 4
if "Hrvatin" in dataset:
    wrap = 3
    
g = sns.FacetGrid(data=components_all,hue="meta",
                  col="imputation",col_wrap=wrap,
                  sharex=False,sharey=False)

g.map(plt.scatter,*colnames,s=size,alpha=.5)
g.add_legend(ncol=2,title="")

for i,ax in enumerate(g.axes.flat):
    ax.set_title(components_all.imputation.unique()[i], fontsize=16, fontweight="bold")
    if i>=wrap:
        ax.set_xlabel("UMAP_1", fontsize=16)
    if i%wrap == 0:
        ax.set_ylabel("UMAP_2", fontsize=16)

for lh in g._legend.legendHandles: 
    lh._sizes = [50] 

plt.savefig("../results/downstream/UMAP_{}.png".format(dataset),dpi=300)
plt.show()



