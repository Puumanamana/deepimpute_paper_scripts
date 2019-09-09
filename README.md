# Scripts to reproduce the figures from DeepImpute paper

DeepImpute is a single cell RNA-seq imputation algorithm. It is available at https://github.com/lanagarmire/deepimpute 

Arisdakessian, Cedric, Olivier Poirion, Breck Yunits, Xun Zhu, and Lana Garmire.
"DeepImpute: an accurate, fast and scalable deep neural network method to impute single-cell RNA-Seq data." bioRxiv (2018): 353607"
https://www.biorxiv.org/content/early/2018/06/22/353607

## Imputing with other methods
The other methods are available in R or python. For the comparison, each imputation results needs to be un-normalized and directly comparable with the initial raw counts.

## Pre-requisites
To use these scripts, you need first to:
- Download the datasets:
  - Jurkat, 293T, neuron9k from the 10X genomics platform
  - GSE67602, GSE99330 (Fish_Dropseq), GSE102827 (Hrvatin dataset)
  - The RNA FISH dataset is avaiable at https://www.dropbox.com/sh/g9c84n2torx7nuk/AABZei_vVpcfTUNL7buAp8z-a?dl=0
- Impute all the datasets with all 7 methods (or the ones that can be ran)
- Undo normalization if the result is normalized
- Collect all results and organize them the following way:
  - Accuracy experiment: All datasets (jurkat, 293T, neuron9k, GSE67602) are stored in a single h5 file, where each dataset is in a different group, and each group keys correspond to all methods + the cell labels (`cells` key), the gene labels (`genes` key), the masked dataset (`raw` key), and the true dataset ('truth' key)
  - FISH experiment: There are 3 groups: `dropseq`, `fish`, `imputed`. the dropseq group consists in the raw dropseq data (`raw`) + the `genes` and `cells` labels. The `fish` group contains the same information in `data`, `genes` and `cells`. The `imputed` group has one key per imputation method.
  - Downstram analysis experiment: One .h5ad (`scanpy` format) per dataset (sim or Hrvatin) per method. The naming is `${method}_${dataset}.h5ad`

# Running the scripts
You can run the scripts inside a Docker container. Docker is available at https://www.docker.com/
Once installed, you need to build the container:
```
# docker build . -t deepimpute_figures
```
Once built, you can access the container using the following command:
```
# docker run -it deepimpute_figures
```

The scripts can then easily be launched using `python script.py` or `Rscript script.R` for R scripts.

## Code
The code is separated into 2 main categories:
- The run_*.py files imputes the data with deepimpute, load the results of other methods (located in `paper_data`), extract some metrics, and save them in the `results` folder
- The scripts in `plot_utils/` simply plot the metrics in the `results` folder.

Some other scripts are available:
- The data masking procedure for the accuracy experiment is in `accuracy/data_masking.py`
- The simulation parameters are available in the `downstream` folder


