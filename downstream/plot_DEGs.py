import pickle
import seaborn as sns
import matplotlib.pyplot as plt

from sklearn.metrics import roc_curve,roc_auc_score
import pandas as pd

import sys
sys.path.append('..')
from config import colors

def load_DEGs():
    with open("../results/downstream/DEGs.pickle","rb") as handle:
        DEGs = pickle.load(handle)
    return DEGs

def load_truth():
    gene_info = pd.read_csv("../paper_data/downstream/gene_info_sim.csv",index_col=0)
        
    nGroups = gene_info.shape[1] - 4
    groups = [ "Group{}".format(g) for g in range(1,nGroups+1) ]
                
    DEGs_true = {group: [ gene for gene in gene_info.index
                          if gene_info.loc[gene,"DEFac{}".format(group)] !=1]
                 for group in groups }
    return DEGs_true

def assess_prediction(DEGs,DEGs_true,method):
    y_hat = []
    y = []

    for group in DEGs:
        probs_adj = 1-DEGs[group]['adj_pval']
        y_hat += probs_adj.tolist()
        y += [ 1 if gene in DEGs_true[group] else 0 for gene in probs_adj.index ]

    data = pd.DataFrame([y_hat,y], index=["y_hat","y"]).T
    data["method"] = method

    return data

def plot_ROC(df):
    ROC_data = {'fpr': [], 'tpr': [], 'method': [] }
    AUCs = {}

    fig,ax = plt.subplots()

    for name, data in df.groupby('method').agg(list).iterrows():
        AUCs[name] = roc_auc_score(data.y, data.y_hat)
        fpr, tpr, thresholds = roc_curve(data.y, data.y_hat, pos_label=1)
        ROC_data['fpr'] += fpr.tolist()
        ROC_data['tpr'] += tpr.tolist()
        ROC_data['method'] += ["{0} AUC={1:.3f}".format(name,AUCs[name])]*len(fpr)

    ROC_data = pd.DataFrame(ROC_data)
    AUCs = pd.Series(AUCs)
    
    sns.lineplot(data=ROC_data,x='fpr',y='tpr',hue='method',drawstyle='steps-pre',ax=ax,
                 palette=colors.drop(['FISH','truth']))
    
    ax.set_title('ROC curves',fontsize=16)
    ax.set_xlabel("FPR",fontsize=16) ; ax.set_ylabel("TPR",fontsize=16)

    handles, labels = ax.get_legend_handles_labels()
    order = AUCs.argsort()
    plt.legend([handles[1:][idx] for idx in order[::-1]],
               [labels[1:][idx] for idx in order[::-1]])
    plt.savefig("../results/downstream/DEGs_sim.png",dpi=300)

    plt.show()
    
if __name__ == '__main__':
    DEGs = load_DEGs()
    DEGs_true = load_truth()
    
    res = [ assess_prediction(df,DEGs_true,method) for method,df in DEGs.items() ]
    res_df = pd.concat(res)

    plot_ROC(res_df)
