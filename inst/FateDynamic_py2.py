from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter
import warnings
warnings.simplefilter("ignore")

import scanpy as sc
import scvelo as scv
import anndata
import palantir
import os
from os import path
import pandas as pd
from scvelo import logging
import numpy as np
import matplotlib
from matplotlib import font_manager
import matplotlib.pyplot as plt
from matplotlib.pyplot import savefig

# Parse command line arguments
parser = ArgumentParser(formatter_class = ArgumentDefaultsHelpFormatter)
parser.add_argument("-i","--input",type=str,help="the address of inputted anndata file")
parser.add_argument("-m","--minCounts",default=20,type=int,help="minimum counts for both spliced and unspliced matrix")
parser.add_argument("-n","--nhvgs",default=2000,type=int,help="the number of highly variable genes")
parser.add_argument("-pc","--pcs",default=30,type=int,help="the number of PCs")
parser.add_argument("-k","--knn",default=30,type=int,help="the number of nearest neighbors")
parser.add_argument("-pl","--palantir",choices = ('True','False'),default='True',type=str,help="using palantir to perform fate prediction")
parser.add_argument("-dc","--dcs",default=10,type=int,help="the number of diffusion components used in palantir")
parser.add_argument("-r","--root",default=None,type=str,help="the id of root cell")
parser.add_argument("-nw","--nwps",default=500,type=int,help="the number of waypoints")
parser.add_argument("-mo","--mode",default="dynamical",type=str,help="the mode of RNA velocity")
parser.add_argument("-p","--plot",choices = ('True','False'),default='True',type=str,help="output the velocity stream embedding plot")
parser.add_argument("-f","--fate",choices = ('True','False'),default='True',type=str,help="if perform cellrank fate prediction")
parser.add_argument("-c","--cluster",default=None,type=str,help="the cell type/cluster identity label")
parser.add_argument("-w","--weight",default=0.2,type=float,help="the weight for building transition matrix in cellrank")
parser.add_argument("-nc","--ncores",default=6,type=int,help="the number of cores")
parser.add_argument("-t","--tolerance",default=1e-6,type=float,help="tolerance")
args = vars(parser.parse_args())

## setup parameters
input = args["input"]
min_counts = args["minCounts"]
nhvgs = args["nhvgs"]
npcs = args["pcs"]
knn = args["knn"]
palantir_model = args["palantir"] == 'True'
ndcs = args["dcs"]
root = args["root"]
num_waypoints = args["nwps"]
mode = args["mode"]
plot = args["plot"] == 'True'
fate_predict = args["fate"] == 'True'
cluster_label = args["cluster"]
weight_connectivities = args["weight"]
ncores = args["ncores"]
tol = args["tolerance"]
logging.info("input:"+str(input))
logging.info("min_counts:"+str(min_counts))
logging.info("nhvgs:"+str(nhvgs))
logging.info("npcs:"+str(npcs))
logging.info("knn:"+str(knn))
logging.info("palantir_model:"+str(palantir_model))
logging.info("ndcs:"+str(ndcs))
logging.info("root:"+root)
logging.info("num_waypoints:"+str(num_waypoints))
logging.info("mode:"+str(mode))
logging.info("plot:"+str(plot))
logging.info("fate_predict:"+str(fate_predict))
logging.info("cluster_label:"+str(cluster_label))
logging.info("weight_connectivities:"+str(weight_connectivities))
logging.info("ncores:"+str(ncores))
logging.info("tol:"+str(tol))

############################ function ################################
def plot_palantir_results(pr_res, tsne,writekey,save,s=3):
    """ Plot Palantir results on tSNE
    """

    # Set up figure
    n_branches = pr_res.branch_probs.shape[1]
    n_cols = 6
    n_rows = int(np.ceil(n_branches / n_cols))
    fig = plt.figure(figsize=[2 * n_cols, 2 * (n_rows + 2)])
    gs = plt.GridSpec(
        n_rows + 2, n_cols, height_ratios=np.append([0.75, 0.75], np.repeat(1, n_rows))
    )
    cmap = matplotlib.cm.plasma
    # Pseudotime
    ax = plt.subplot(gs[0:2, 1:3])
    c = pr_res.pseudotime[tsne.index]
    ax.scatter(tsne.iloc[:, 0], tsne.iloc[:, 1], s=s, cmap=matplotlib.cm.plasma, c=c)
    normalize = matplotlib.colors.Normalize(vmin=np.min(c), vmax=np.max(c))
    cax, _ = matplotlib.colorbar.make_axes(ax)
    cbar = matplotlib.colorbar.ColorbarBase(cax, norm=normalize, cmap=cmap)
    ax.set_axis_off()
    ax.set_title("Pseudotime")

    # Entropy
    ax = plt.subplot(gs[0:2, 3:5])
    c = pr_res.entropy[tsne.index]
    ax.scatter(tsne.iloc[:, 0], tsne.iloc[:, 1], s=s, cmap=matplotlib.cm.plasma, c=c)
    normalize = matplotlib.colors.Normalize(vmin=np.min(c), vmax=np.max(c))
    cax, _ = matplotlib.colorbar.make_axes(ax)
    cbar = matplotlib.colorbar.ColorbarBase(cax, norm=normalize, cmap=cmap)
    ax.set_axis_off()
    ax.set_title("Differentiation potential")

    for i, branch in enumerate(pr_res.branch_probs.columns):
        row = int(np.floor(i / n_cols))
        ax = plt.subplot(gs[row + 2, np.remainder(i, n_cols)])
        c = pr_res.branch_probs.loc[tsne.index, branch]
        ax.scatter(
            tsne.iloc[:, 0], tsne.iloc[:, 1], s=s, cmap=matplotlib.cm.plasma, c=c
        )
        normalize = matplotlib.colors.Normalize(vmin=np.min(c), vmax=np.max(c))
        cax, _ = matplotlib.colorbar.make_axes(ax)
        cbar = matplotlib.colorbar.ColorbarBase(cax, norm=normalize, cmap=cmap)
        ax.set_axis_off()
        ax.set_title(branch, fontsize=10)
    if save:
        savefig(writekey, dpi=300)
#############################################################################

if __name__ == "__main__":
    if palantir_model:
        logging.info("Using palantir to perform cell fate prediction")
        ad = sc.read_h5ad(input)
        sc.pp.filter_genes(ad, min_counts=min_counts)
        sc.pp.normalize_per_cell(ad)
        sc.pp.log1p(ad)
        sc.pp.highly_variable_genes(ad, n_top_genes=nhvgs)
        sc.pp.pca(ad)
        pca_projections = pd.DataFrame(ad.obsm['X_pca'], index=ad.obs_names)
        dm_res = palantir.utils.run_diffusion_maps(pca_projections, n_components=ndcs)
        ms_data = palantir.utils.determine_multiscale_space(dm_res)
        pr_res = palantir.core.run_palantir(ms_data, root, num_waypoints=num_waypoints)
        pr_res.branch_probs.columns = ad.obs[cluster_label][pr_res.branch_probs.columns].values
        
        fate_prob = pr_res._branch_probs
        pseudotime = pr_res._pseudotime
        if plot:
            sc.pp.neighbors(ad, n_pcs=npcs, n_neighbors=knn)
            sc.tl.umap(ad)
            umap = pd.DataFrame(ad.obsm['X_umap'], index=ad.obs_names,columns = ["x","y"])
            if path.exists("palantir_figures/"):
                plot_palantir_results(pr_res = pr_res, tsne = umap,writekey="palantir_figures/palantir_res.png",save=True)
                sc.pl.embedding(ad, basis='umap',color = cluster_label,legend_loc='on data', title='',legend_fontsize='x-small',save='.png')
            else:
                os.mkdir("palantir_figures/")
                plot_palantir_results(pr_res = pr_res, tsne = umap,writekey="palantir_figures/palantir_res.png",save=True)
                sc.pl.embedding(ad, basis='umap',color = cluster_label,legend_loc='on data', title='',legend_fontsize='x-small',save='.png')
                
        logging.info("Identify " + str(fate_prob.shape[1])+" terminal states")
        logging.info("saving...")
        sampleID = ad.obs.index.tolist()
        geneID = ad.var.index.tolist()
        
        if isinstance(ad.X,np.ndarray):
            raw = pd.DataFrame(ad.X.T,index=geneID,columns=sampleID)
        else:
            raw = pd.DataFrame(ad.X.toarray().T,index=geneID,columns=sampleID)
        
        adata_save = anndata.AnnData(X = raw.T,obs = fate_prob)
        adata_save.obs["pseudotime"] = pseudotime
    else:
        adata = sc.read_h5ad(input)
        scv.pp.filter_genes(adata, min_shared_counts=min_counts)
        scv.pp.normalize_per_cell(adata)
        adata_raw = adata.copy()
        scv.pp.filter_genes_dispersion(adata, n_top_genes=nhvgs)
        scv.pp.log1p(adata)
        scv.pp.log1p(adata_raw)

        sc.tl.pca(adata)  
        sc.pp.neighbors(adata, n_pcs=npcs, n_neighbors=knn)
        sc.tl.umap(adata)
        scv.pp.moments(adata,n_pcs=npcs,n_neighbors=knn)
        ## assign the pca and neighborhood graph
        adata_raw.obsm['X_pca'] = adata.obsm['X_pca'].copy()
        adata_raw.obsp['distances'] = adata.obsp['distances'].copy()
        adata_raw.obsp['connectivities'] = adata.obsp['connectivities'].copy()
        adata_raw.uns['pca'] = adata.uns['pca'].copy()
        adata_raw.uns['neighbors'] = adata.uns['neighbors'].copy()
        scv.pp.moments(adata_raw,n_pcs=npcs,n_neighbors=knn)

        logging.info("Using " + mode + " mode to perform velocity analysis...")
        if mode == "dynamical":
            logging.info("Estimating RNA velocity for HVGs...")
            scv.tl.recover_dynamics(adata,n_jobs = ncores)
            scv.tl.velocity(adata,mode='dynamical',n_jobs = ncores)
            logging.info("Estimating RNA velocity for All Genes...")
            scv.tl.recover_dynamics(adata_raw,n_jobs = ncores)
            scv.tl.velocity(adata_raw,mode='dynamical',n_jobs = ncores)
        else:
            logging.info("Estimating RNA velocity for HVGs...")
            scv.tl.velocity(adata,mode=mode,n_jobs = ncores)
            logging.info("Estimating RNA velocity for All Genes...")
            scv.tl.velocity(adata_raw,mode=mode,n_jobs = ncores)

        scv.tl.velocity_graph(adata)
        scv.tl.velocity_pseudotime(adata)

        if plot:
            scv.pl.velocity_embedding_stream(adata, color=cluster_label,figsize=[6,4],dpi=300,save="velocity_embedding_stream.png",legend_loc='right margin')
                
        if fate_predict:
            import cellrank as cr
            k = cr.tl.transition_matrix(
                adata, weight_connectivities=weight_connectivities, softmax_scale=4, show_progress_bar=False
            )
            g = cr.tl.estimators.GPCCA(k)
            g.compute_schur()
            g.compute_macrostates(cluster_key=cluster_label)
            g.compute_terminal_states(method = "eigengap")
            g.compute_absorption_probabilities(n_jobs = ncores)
            fate_prob = g.absorption_probabilities
            sampleID = adata.obs.index.tolist()
            fate_name = fate_prob.names.tolist()
            fate_prob = pd.DataFrame(fate_prob,index= sampleID,columns=fate_name)

        logging.info("Identify " + str(fate_prob.shape[1])+" terminal states")
        logging.info("saving...")
        sampleID = adata_raw.obs.index.tolist()
        geneID = adata_raw.var.index.tolist()
        veloGene = adata_raw.var["velocity_genes"].tolist()
        #raw = pd.DataFrame(adata_raw.layers["spliced"].T,index=geneID,columns=sampleID)
        
        if isinstance(adata_raw.X,np.ndarray):
            raw = pd.DataFrame(adata_raw.X.T,index=geneID,columns=sampleID)
        else:
            raw = pd.DataFrame(adata_raw.X.toarray().T,index=geneID,columns=sampleID)
            
        Ms = pd.DataFrame(adata_raw.layers["Ms"].T,index=geneID,columns=sampleID)
        velocity = pd.DataFrame(adata_raw.layers["velocity"].T,index=geneID,columns=sampleID)
        #raw = raw[veloGene]
        #Ms = Ms[veloGene]
        #velocity = velocity[veloGene]

        ### Build adata object and save
        var = pd.DataFrame(veloGene,index=geneID,columns=["velocity_genes"])
        adata_save = anndata.AnnData(X = raw.T,obs = fate_prob)
        adata_save.obs["velocity_pseudotime"] = adata.obs["velocity_pseudotime"]
        adata_save.layers["Ms"] = Ms.T
        adata_save.layers["velocity"] = velocity.T
        adata_save.var = var
        
    if palantir_model:
        if path.exists("output/"):
            adata_save.write("output/FateRes.h5ad")
        else:
            os.mkdir("output/")
            adata_save.write("output/FateRes.h5ad")
    else:
        if path.exists("output_velo/"):
            adata_save.write("output_velo/FateRes.h5ad")
        else:
            os.mkdir("output_velo/")
            adata_save.write("output_velo/FateRes.h5ad")
