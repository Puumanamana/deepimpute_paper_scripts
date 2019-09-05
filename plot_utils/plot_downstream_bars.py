import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import sys

from config import colors

dataset = "sim"
clustering_method = "leiden"
METRICS = [ "adjusted_rand_score",
            "adjusted_mutual_info_score",
            "Fowlkes-Mallows",
            "silhouette_score"]

if len(sys.argv)>1:
    dataset = sys.argv[1]

filename = "../results/downstream/{}_{}_clustering_scores.csv".format(clustering_method,dataset)

scores = pd.read_csv(filename,index_col=0)

scores = scores[METRICS]
scores.columns = [score.replace("_","\n") for score in scores.columns]

scores = pd.melt(scores.reset_index(),
                 id_vars="imputation")

methods = [met for met in colors.index if met in np.unique(scores['imputation'])]

sns.barplot(x="variable",y="value",hue="imputation",data=scores,
            palette=colors.loc[methods],
            hue_order=methods)
plt.xlabel("")
plt.ylabel("Score", fontsize=10)
plt.xticks(fontsize=10)

plt.legend(fontsize=10)
plt.title(dataset, fontsize=10, fontweight="bold")

plt.savefig("../results/downstream/barplot_{}.png".format(dataset),dpi=300)
plt.show()
