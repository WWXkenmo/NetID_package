# NetID_package
A scalable method to infer fate specific gene regulatory network from single cell gene expression data

<img src="https://github.com/WWXkenmo/NetID_package/blob/figures/figures/Concept_fig1.png" alt="NetID" width="600" />

## Installation
### Create conda environment
```
conda create --name NetID python=3.8 r-essentials r-base
conda activate NetID # if it encounter the error, run 'source activate' ahead of this code
```
To speed up installation, user could use conda install mamba at first, then use mamba to install other modules
```
conda install mamba
mamba install -c bioconda -c conda-forge cellrank-krylov
mamba install -c r r-devtools
mamba install -c r r-fateid
```
finally, install palantir package
```
pip install palantir
```

### Install RaceID package
clone the git repo into local directory
```
git clone https://github.com/WWXKenmo/NetID_package.git
cd NetID_package
```
install RaceID dependent packages
```
install.packages(c('coop', 'FNN', 'fpc', 'ggplot2', 'harmony', 'ica', 'igraph', 'irlba', 'leiden', 'locfit', 'matrixStats', 'pheatmap', 'princurve', 'quadprog', 'randomForest', 'runner', 'RColorBrewer', 'Rtsne', 'umap', 'vegan'))
```
Then install RaceID in R
```
install.packages("RaceID_0.2.9.tar.gz",repos=NULL, type="source",INSTALL_opts=c("--no-multiarch"))
```

### Install NetID package in R
```
install_github("WWXKenmo/NetID_package")
```
or
```
install.packages("NetID_0.1.0.tar.gz",repos=NULL, type="source",INSTALL_opts=c("--no-multiarch"))
```
