import pandas as pd
import numpy as np
from scipy.stats import pearsonr,ks_2samp

import seaborn as sns
import matplotlib.pyplot as plt

import os,sys
PARENT_DIR = os.path.join(sys.path[0], '..')
sys.path.insert(1, PARENT_DIR)

output_dir = "{}/results/fish".format(PARENT_DIR)

from config import colors

def plot_gini(df):
    color_orders = colors.drop(["DrImpute","truth","FISH"])
    df_melt = pd.melt(df,id_vars=["FISH","Gene"])
    
    g = sns.FacetGrid(df_melt,
                      col="variable",col_wrap=4, col_order=color_orders.index,
                      hue="Gene",
                      palette=color_orders)
    g.map(sns.scatterplot, "FISH", "value",edgecolor="k", s=100, alpha=0.8).add_legend()

    plt.setp(g._legend.get_title(), fontweight='bold')
    scores = {'Pearson': {},'MSE':{}}
    
    for i,ax in enumerate(g.axes.flat):
        ax.plot((0,1), (0,1), c=".2", ls="--")

        ax.set_xlabel('Gini coefficient (FISH)',fontsize=12)
        method = color_orders.index[i]

        if method.lower() == 'raw':
            ax.set_ylabel('Gini coefficient (Raw)',fontsize=12)
        else:
            ax.set_ylabel('Gini coefficient (Imputed)',fontsize=12)

        scores['Pearson'][method] = pearsonr(df['FISH'],df[method])[0]
        scores['MSE'][method] = np.mean((np.array(df['FISH'])-np.array(df[method]))**2)

        ax.text(0.4,0.1,r"$MSE={0:.3g}$".format(scores['MSE'][method]) + "\n" +
                        r"$correlation={0:.3f}$".format(scores['Pearson'][method]),
                fontsize=10)
        ax.set_title(ax.get_title().replace("variable = ",""), fontweight='bold')

    plt.savefig("{}/gini_scatter.png".format(output_dir),dpi=300)
    plt.show()


def plot_distributions(distributions):
    color_orders = colors.drop(["DrImpute","truth"])
    gene_names = distributions.gene.unique()

    g = sns.FacetGrid(distributions, hue="method", hue_order=color_orders.index,
                      col="gene", col_wrap=3,
                      sharey=False, sharex=False,
                      col_order=gene_names)

    g.map(sns.kdeplot, "value").add_legend()
    plt.setp(g._legend.get_title(), fontweight='bold')

    for i, ax in enumerate(g.axes.flat):
        fish_gene = distributions.value[ (distributions.method=="FISH") & (distributions.gene==gene_names[i]) ]
        ax.set_xlim([0, np.percentile(fish_gene,90)])

        if gene_names[i] == "VGF":
            ax.set_ylim([0, 0.3])
        elif gene_names[i] == "TXNRD1":
            ax.set_ylim([0, 0.1])
        
        ax.set_title(gene_names[i],fontsize=16,fontweight='bold')

        if i%3 == 0:
            ax.set_ylabel("Density", fontsize=14)
        if int(i/3) == 1 or i==2:
            ax.set_xlabel("Normalized values", fontsize=14)
    plt.savefig("{}/distributions_FISH.png".format(output_dir),dpi=300)
    plt.show()

def get_ks_values(data_distr):
    grouped = data_distr.groupby(['method','gene']).agg(list)
    results = []
    
    for method in data_distr.method.unique():
        if method != 'FISH':
            for gene in data_distr.gene.unique():
                imputed = grouped.loc[(method,gene),'value']
                fish = grouped.loc[('FISH',gene),'value']
                results.append([ks_2samp(fish,imputed).statistic,
                                method,gene])

    results = pd.DataFrame(results,columns=['KS_score','method','gene'])
    results = results.pivot('method','gene')
    results = results.iloc[np.argsort(results.sum(axis=1)).values,:]
    
    results.to_csv('{}/KS_stats.csv'.format(output_dir))
        
if __name__ == '__main__':
    data_gini = pd.read_csv("{}/gini.csv".format(output_dir)).dropna()
    plot_gini(data_gini)
    
    data_distr = pd.read_csv("{}/normalized_distributions.csv".format(output_dir), index_col=0)
    plot_distributions(data_distr)
    get_ks_values(data_distr)
