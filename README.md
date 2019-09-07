# Scripts to reproduce the figures from DeepImpute paper

DeepImpute is a single cell RNA-seq imputation algorithm. It is available at https://github.com/lanagarmire/deepimpute 

Arisdakessian, Cedric, Olivier Poirion, Breck Yunits, Xun Zhu, and Lana Garmire.
"DeepImpute: an accurate, fast and scalable deep neural network method to impute single-cell RNA-Seq data." bioRxiv (2018): 353607"
https://www.biorxiv.org/content/early/2018/06/22/353607

To use these scripts, you need to first impute the paper's datasets using the methods that are compared (SAVER, scImpute, MAGIC, DCA, DrImpute and VIPER). The matrices are then stored in a h5 file for speed and clarity purposes. However, for the downstream analysis figure, keeping the files on the h5ad format (scanpy) facilitates the analysis.

# Running the scripts
You can run the scripts inside a Docker container. Docker is available at https://www.docker.com/
Once installed, you need to build the container:
```
# docker build . -t deepimpute_figures
```
To run the scripts inside the container, you just need to do
```
# docker run -it deepimpute_figures
```

# Code
The code is separated into 2 main categories:
- The run_*.py files imputes the data with deepimpute, load the results of other methods (located in `paper_data`), extract some metrics, and save them in the `results` folder
- The scripts in `plot_utils/` simply plot the metrics in the `results` folder.

Some other scripts are available:
- The data masking procedure for the accuracy experiment is in `accuracy/data_masking.py`
- The simulation parameters are available in the `downstream` folder
