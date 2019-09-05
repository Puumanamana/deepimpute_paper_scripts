import pandas as pd
import seaborn as sns
from matplotlib import rcParams

PLOT_PARAMS = {"title":{"fontweight": "bold","fontsize": 18},
               "axes":{"fontweight": "bold","fontsize": 18},
               "heatmap":{"aspect": "auto", "vmax": 20,"cmap": "Greens"}}

rcParams['font.family'] = "serif"

imputation_methods = ["deepImpute","DCA","MAGIC","SAVER","scImpute","VIPER","DrImpute"]
datasets_order = [ "jurkat", "293T", "neuron9k", "GSE67602" ]

plot_order = imputation_methods + ["raw", "FISH", "truth"]

colors = pd.Series({method: sns.color_palette()[i]
                    for i,method in enumerate(plot_order)
                    })
