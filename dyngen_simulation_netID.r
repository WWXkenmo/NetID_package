library(tidyverse)
library(dyngen)

######## Figure 3 simulation ############## 
backbone <- backbone_bifurcating()
config <- 
    initialise_model(
        backbone = backbone,
        num_cells = 4000,
        num_tfs = 50,
        download_cache_dir = tools::R_user_dir("dyngen","data"),
        num_targets = 200,
        num_hks = 50,
        verbose = FALSE
    )

out <- generate_dataset(
    config,
    format = "anndata",
    make_plots = TRUE
)

######## Figure 3 simulation ############## 

backbone <- backbone_cycle()
config <- 
    initialise_model(
        backbone = backbone,
        num_cells = 4000,
        num_tfs = 50,
        download_cache_dir = tools::R_user_dir("dyngen","data"),
        num_targets = 200,
        num_hks = 50,
        verbose = FALSE
    )

out <- generate_dataset(
    config,
    format = "anndata",
    make_plots = TRUE
)

######## Figure 5 simulation ############## 
backbone <- backbone_bifurcating()
config <- 
    initialise_model(
        backbone = backbone,
        num_cells = 4000,
        num_tfs = 50,
        simulation_params = simlation_default(compute_cellwise_grn = TRUE),
        num_targets = 200,
        num_hks = 50,
        verbose = FALSE
    )

out <- generate_dataset(
    config,
    format = "anndata",
    make_plots = TRUE
)
