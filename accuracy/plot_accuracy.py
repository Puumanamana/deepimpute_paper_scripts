#!/usr/bin/env python

import matplotlib.pyplot as plt
from matplotlib.colors import Normalize

import numpy as np
import pandas as pd
import seaborn as sns
from os.path import dirname

from sklearn.metrics import mean_squared_error
from scipy.stats import pearsonr

import sys
sys.path.append('..')
from config import colors, datasets_order

COLORS = colors.drop(["raw","truth","FISH"])
METHODS = COLORS.index

dims = ["cell level","gene level"]
outputdir = "../results/accuracy"

def mse_boxplot(folder=outputdir):

    files = [ "{}/MSE_{}.csv".format(folder,dataset) for dataset in datasets_order ]
    data = pd.concat([pd.read_csv(f,index_col=0) for f in files])
    data.method = [ x.upper() if x=='dca' else x for x in data.method ]

    fig = sns.catplot(x="dataset",y="MSE",hue="method",hue_order=METHODS,row="axis",row_order=['cells','genes'],
                      kind="box",data=data,showfliers=False,legend=False,palette=COLORS)

    for i,ax in enumerate(fig.axes.flat):
        ax.set_title("")
        ax.set_xlabel("")
        ax.set_ylabel("MSE ({})".format(dims[i]), fontsize=14, fontweight="bold")
        plt.setp(ax.get_xticklabels(), fontsize=14, fontweight="bold", rotation=15)
    
    plt.legend(loc="best")

    plt.savefig("{}/boxplot.png".format(folder),dpi=300)

def scatter_plot(folder=outputdir,lims=[0,2.5]):

    files = [ "{}/scatter_{}.csv".format(folder,dataset) for dataset in datasets_order ]    
    data = pd.concat([pd.read_csv(f,index_col=0) for f in files])
    data.method = [ x.upper() if x=='dca' else x for x in data.method ]    

    g = sns.FacetGrid(data, row="dataset", col="method", despine=True)
    g.map(sns.scatterplot,"Truth","Imputed",s=20,alpha=.1)

    for i, ax in enumerate(g.axes.flat):
        ax.plot(lims,lims,'r-.')
        ax.set_xlim(lims) ; ax.set_ylim(lims)
        ax.set_xlabel("") ; ax.set_ylabel("")
        ax.set_title("", fontsize=10)

        method = METHODS[i%len(METHODS)]
        dataset = datasets_order[int(i/len(METHODS))]
        data_i = data[(data['dataset']==dataset) & (data['method']==method)]

        if len(data_i) >0 :
            mseScore = mean_squared_error(data_i["Truth"].values,
                                          data_i["Imputed"].values)
            pearsonScore = pearsonr(data_i["Truth"].values,
                                    data_i["Imputed"].values)[0]
            
            title = "MSE = {0:.3f}\nPearson = {1:.3f} \n".format(mseScore, pearsonScore)
            ax.text(0.5,1.8,title)
        
        if i < len(METHODS):
            ax.set_title(METHODS[i], fontsize=16,fontweight='bold')
        if i % len(METHODS) == 0:
            ax.set_ylabel(datasets_order[int(i/len(METHODS))], fontsize=16,fontweight='bold')
            
    plt.savefig("{}/scatterplot.png".format(folder),dpi=300)

def display_table(df):
    fig,ax = plt.subplots()
    ax.axis('off')

    table = ax.table(
        cellText=np.array([["{:.2f}".format(col)
                            for col in row]
                           for row in df.values]),
        rowLabels=df.index,
        colLabels=df.columns,
        colWidths = [0.2]*df.shape[0],
        loc='center',
        cellColours=plt.cm.Greens(Normalize(1,300)(df.values)) )
    table.auto_set_font_size(False)
    table.set_fontsize(10)

def pvalue_plot(filename="{}/pvalues_table.csv".format(outputdir)):

    p_values = pd.read_csv(filename,index_col=0)
    p_values.method = [ x.upper() if x=='dca' else x for x in p_values.method ]
    
    for axis in ["cell","gene"]:
        fig,ax = plt.subplots()
        
        data = p_values[p_values.axis==axis].drop("axis",axis=1)
        data = -np.log10(data.pivot("dataset","method"))
        data.columns = data.columns.levels[1]
        data = data.reindex(index=datasets_order,columns=METHODS[1:])

        display_table(data.T)
        plt.savefig("{}/pvalues_{}.png".format(dirname(filename),axis),dpi=300)        
            
if __name__ == '__main__':

    mse_boxplot()
    scatter_plot()
    pvalue_plot()
    plt.show()


