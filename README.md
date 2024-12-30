# NetID
A scalable method to infer fate specific gene regulatory network from single cell gene expression data

<img src="https://github.com/WWXkenmo/NetID_package/blob/figures/figures/Concept_fig1.png" alt="NetID" width="600" />

## Tutorial

[Tutorial of NetID](https://htmlpreview.github.io/?https://github.com/WWXkenmo/NetID_package/blob/main/inst/doc/NetID.html)

## Installation
### Basic installation
#### Create conda environment (recommand but not necessary)
```
conda create --name NetID -c conda-forge -c bioconda r-seurat=4* python=3.10 r-essentials r-base=4.2.0
conda activate NetID # if it encounter the error, run 'source activate' ahead of this code
```
#### Install devtools and geosketch
```
conda install -c conda-forge r-devtools
pip install geosketch
```
#### install dependence in R including GENIE3, SingleCellExperiment and NetID (switch to the R terminal)
```
install.packages("BiocManager")
BiocManager::install("GENIE3")
BiocManager::install("SingleCellExperiment")
devtools::install_github("WWXKenmo/NetID_package")
```
### Advance installation
#### install cellrank and palantir to realize lineage-specific GRN prediction
To speed up installation, user could use conda install mamba at first, then use mamba to install other modules
```
conda install mamba -c conda-forge
mamba install -c bioconda -c conda-forge cellrank-krylov
## or could use conda to install
## conda install -c bioconda -c conda-forge cellrank-krylov

## install proper version of package
pip install scanpy==1.9.2
pip install matplotlib==3.7
pip install pandas==1.5.3
pip install palantir==1.0.1
pip uninstall numpy
pip install numpy==1.23.5
```
#### install cytotrace and scent to determine the root cell
```
devtools::install_github("aet21/SCENT")
```

### Citation
```
@article{PUSH:72058,
  author = {Wang, W. and Wang, Y. and Lyu, R. and Gr{\"u}n, D.},
  title = {{Scalable identification of lineage-specific gene regulatory networks from metacells with NetID.}},
  journal = {Genome Biol.},
  location = {Campus, 4 Crinan St, London N1 9xw, England},
  publisher = {Bmc},
  volume = {25},
  number = {1},
  year = {2024},
  issn = {1474-760X},
  eissn = {1465-6906},
}
```
