# NetID_package
A scalable method to infer fate specific gene regulatory network from single cell gene expression data

<img src="https://github.com/WWXkenmo/NetID_package/blob/figures/figures/Concept_fig1.png" alt="NetID" width="600" />

## Tutorial

```
vignette("NetID")
```

## Installation
### Basic installation
#### Create conda environment (recommand but not necessary)
```
conda create --name NetID python=3.10 r-essentials r-base=4.2.0
conda activate NetID # if it encounter the error, run 'source activate' ahead of this code
```
#### Install devtools and geosketch
```
conda install -c conda-forge r-devtools
pip install geosketch
```
#### install NetID
```
install.packages('NetID_0.1.0.tar.gz', repos=NULL, type='source')
```
### Advance installation
#### install cellrank and palantir to realize lineage-specific GRN prediction
To speed up installation, user could use conda install mamba at first, then use mamba to install other modules
```
conda install mamba -c conda-forge
mamba install -c bioconda -c conda-forge cellrank-krylov
## or could use conda to install
## conda install -c bioconda -c conda-forge cellrank-krylov
pip install numpy==1.23.5 palantir
```
#### install summa to output global GRN
```
devtools::install_github("learn-ensemble/R-SUMMA")
```
#### install cytotrace and scent to determine the root cell
```
devtools::install_github("aet21/SCENT")
