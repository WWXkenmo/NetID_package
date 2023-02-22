# NetID_package
A scalable method to infer fate specific gene regulatory network from single cell gene expression data

<img src="https://github.com/WWXkenmo/NetID_package/blob/figures/figures/Concept_fig1.png" alt="NetID" width="600" />

## Installation
### Create conda environment
```
conda create --name NetID python=3.8
conda activate NetID # if it encounter the error, run 'source activate' ahead of this code
conda install -c bioconda -c conda-forge cellrank-krylov
pip install palantir
```
To speed up installation, user could use conda install mamba at first, then use mamba to install other modules
```
conda install mamba
mamba install -c bioconda -c conda-forge cellrank-krylov
```

### Install RaceID package
```
install.packages("RaceID_0.2.9.tar.gz",repos=NULL, type="source",INSTALL_opts=c("--no-multiarch"))
```

### Install NetID package in R
```
install_github("WWXKenmo/NetID_package")
```
