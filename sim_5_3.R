# setwd("/Users/gregfaletto/Documents/GitHub/css-sims")
# setwd("/Users/gregfaletto/Google Drive/Data Science/LaTeX/Generalized Stability Selection Presentation")


# plot_eval(sim, "sc4linear")
# plot_evals(sim, "size", "sc4linear")
# sim <- load_simulation("sparse-blocked-linear-model")

# subset_simulation(toy_model, methods=ev_6e@method_name[1:6]) %>%
#     plot_evals("size", "sc4linear", method_col=(c(1, rep(2, 5)))) +
#     scale_color_manual(values=c("red", rep("blue", 5)))

# subset_simulation(toy_model, methods=ev_6e@method_name[1:6]) %>%
#     plot_evals("size", "sc4linear", method_col=(c(1, rep(2, 5)))) + 
#     scale_color_manual(values=c("red", rep("blue", 5)))

# 
# subset_simulation(toy_model, methods=ev_6e@method_name[1:6]) %>%
#      plot_evals("size", "sc4linear") +
#      scale_color_manual(values=c("red", rep("blue", 5)))+geom_path()

# rm(list=ls()[!ls() %in% c("sim_6e")])
# dev.off()
rm(list=ls())


############# Load Libraries #################

library(cssr)

library(simulator) # this file was created under simulator version 0.2.0
library(MASS)
library(glmnet)
library(digest)
library(knitr)
library(ggplot2)
library(doParallel)
library(ccaPP)
library(Metrics)
library(gridExtra)
library(irlba)
library(doMC)
library(parallel)
library(dplyr)
library(cowplot)


# PMA requires library "impute" which is no longer on CRAN. Instructions
# to download as listed on 
# http://www.bioconductor.org/packages/release/bioc/html/impute.html:

# ## try http:// if https:// URLs are not supported
# source("https://bioconductor.org/biocLite.R")
# biocLite("impute")

library(PMA)
library(e1071)

############# Load Functions #################

wd <- getwd()
code_dir <- paste(wd, "Helper Functions", sep="/")
# setwd("/Users/gregfaletto/Google Drive/Data Science/R/USC/Stability Selection/toy_example/New method")
# source(file="stabsel_copy.R")
# wd <- "/Users/gregfaletto/Google Drive/Data Science/R/USC/Stability Selection/toy_example"
# setwd(wd)


registerDoParallel()
registerDoMC(cores = detectCores())

setwd(code_dir)
source("toy_ex_slide_funcs.R")
source("model_functions.R")
source("method_functions.R")
source("eval_functions.R")
setwd(wd)



###############################################

# Run simulation, or load simulation that has been previously run?
run_new_sim <- TRUE
# Titles on p_hat plots?
p_hat_titles <- TRUE
# legends on mse/false selection/proxy selection plots
mse_legend <- TRUE
# Only print plot for specific values? (NA if not)
# p_plot <- 100
p_plot <- NA
# beta_plot <- 2
beta_plot <- NA
# Method (choose MB or SS)
method <- "SS"
if(!(method %in% c("SS", "MB"))){
    stop("!(method %in% c(SS, MB))")
}
# Seeds (first for model simulations, second for metrics)
seed1 <- 457335
seed2 <- 734355
n_model <- 200
# Number of unlabeled observations used for estimating clusters
n_clus <- 200
# Cutoff for absolute correlation for estimated clusters
est_clus_cutoff <- 0.5
n_test <- 10000
n_sims <- 2000
# n_sims <- 10
# p <- 50
p <- 100
nblocks <- 1
sig_blocks <- 1
k_unblocked <- 10
beta_low <- 1
# beta_high <- 2
beta_high <- 1.5
block_size <- 15
n_strong_block_vars <- 5
rho_high <- 0.9
rho_low <- 0.5 
var <- 1
snr <- 3
sigma_eps_sq <- NA
# Number of parameter combinations
n_param_combs <- length(p)*length(n_model)*length(k_unblocked)*length(beta_high)
# Tau for evaluations after simulations
tau <- 0.51
# Maximum number of features to include in stability selection models
p_max <- k_unblocked + sig_blocks
# Methods to evaluate in visualizations
if(method=="SS"){
    methods_to_eval <- c(
        # "subspace_sel_grass_random", 
        # "subspace_ss_max_cancor_random", 
        # "subspace_ss_min_cancor_random",
        # "SS_SS_random",
        # "SS_SS_random_custom",
        # "SS_GSS_random"
        "SS_GSS_random_custom"
        # , "lasso_random"
        # , "SS_GSS_random_avg"
        , "SS_GSS_random_avg_custom"
        # , "SS_GSS_random_avg_unwt"
        , "SS_GSS_random_avg_unwt_custom"
        , "BRVZ_avg_unwt"
        # , "lasso_proto"
        )
}else if(method=="MB"){
    methods_to_eval <- c("lasso_random", "lassoMB_phat_random")
}
n_meths_to_eval <- length(methods_to_eval)

# Stability selection-type methods
stab_meths_to_eval <- intersect(methods_to_eval, c("lassoMB_phat_random",
    "lassoMB_phat_ideal_random", "lassoMB_phat_cor_squared_random",
    "subspace_ss_max_cancor_random", "subspace_sel_grass_random",
    "SS_SS_random", "SS_SS_random_custom", "SS_GSS_random",
    "SS_GSS_random_custom", "lassoSS_phat_cor_squared_random",
    "SS_GSS_random_avg", "SS_GSS_random_avg_custom", "SS_GSS_random_avg_unwt",
    "SS_GSS_random_avg_unwt_custom"))

n_stab_meths_to_eval <- length(stab_meths_to_eval)

non_stab_meths <- setdiff(methods_to_eval, stab_meths_to_eval)

# # Display names of possible methods
# DISPLAY_NAMES <- ("Lasso", "Stability Selection", "Sparse CSS",
#     "Weighted Averaged CSS", "Simple Averaged CSS",
#     "Cluster Representative Lasso", "Protolasso")

# if(!("lasso_random" %in% methods_to_eval)){
#     stop("!(lasso_random %in% methods_to_eval)")
# }

known_cluster_meths <- c("SS_CSS_sparse_cssr" # Sparse cluster stability selection
    # , SS_GSS_random_custom # Sparse cluster stability selection
    , "SS_CSS_weighted_cssr" # Weighted averaged cluster stability
    # selection
    , "SS_CSS_avg_cssr" # Simple averaged cluster stability
    # selection
    , "clusRepLasso_cssr" # Cluster representative lasso
    , "protolasso_cssr" # Protolasso
    )

est_cluster_meths <- c("SS_CSS_sparse_cssr_est" # Sparse cluster stability selection,
    # estimated clusters
    , "SS_CSS_weighted_cssr_est" # Weighted averaged cluster stability
    # selection, estimated clusters
    , "SS_CSS_avg_cssr_est" # Simple averaged cluster stability
    # selection, estimated clusters
    , "clusRepLasso_cssr_est" # Cluster representative lasso, estimated
    # clusters
    , "protolasso_cssr_est"
    )

# Calculate clustered stability metric?

calc_clust_stab <- TRUE

if(calc_clust_stab & (sig_blocks != 1)){
    stop("sig_blocks != 1 (can't calculate clustered stability metric)")
}

setwd(wd)

t0 <- Sys.time()

if(run_new_sim){

    print("Starting simulations...")

    set.seed(seed1)

    if(method=="MB"){
        gss_random_weighted_custom <- new_simulation("gss_random_weighted_custom",
        "GSS (Random Design, Weighted Averaging)") %>%
        generate_model(make_model=make_blocked_lin_mod4_ran_weight,
        n = n_model,
        n_test = n_test,
        p = p,
        # p = as.list(c(2*n_model, 40*n_model)),
        k_unblocked = k_unblocked,
        est_clus_cutoff=est_clus_cutoff,
        beta_low = beta_low,
        beta_high = beta_high,
        nblocks = nblocks,
        sig_blocks = sig_blocks,
        block_size = block_size,
        n_strong_block_vars = n_strong_block_vars,
        rho_high = rho_high,
        rho_low = rho_low,
        var = var,
        snr = snr 
        # , vary_along = "p"
        ) %>% simulate_from_model(nsim = n_sims) %>%
        run_method(c(lassoMB_phat_random
            , lassoMB_phat_ideal_random
            , lassoMB_phat_cor_squared_random
            , lasso_random
            # , subspace_ss_max_cancor_random
            # , list_of_lasso_ps_sss_5_.7_random
            # , list_of_lasso_sss_5_.51_random
            # , lassoSS_random
        )) %>% evaluate(list(phat, labels))
    } else{
        gss_random_weighted_custom <- new_simulation("gss_random_weighted_custom",
        "GSS (Random Design, Weighted Averaging)") %>%
        generate_model(make_model=make_blocked_lin_mod4_ran_weight,
        n = n_model,
        n_clus = n_clus,
        n_test = n_test,
        p = p,
        k_unblocked = k_unblocked,
        est_clus_cutoff=est_clus_cutoff,
        beta_low = beta_low,
        beta_high = beta_high,
        nblocks = nblocks,
        sig_blocks = sig_blocks,
        block_size = block_size,
        n_strong_block_vars = n_strong_block_vars,
        rho_high = rho_high,
        rho_low = rho_low,
        var = var,
        snr = snr 
        # , vary_along = ifelse(n_param_combs > 1, "p", NULL)
        # , vary_along = c("p", "beta_high")
        ) %>% simulate_from_model(nsim = n_sims)

        save_simulation(gss_random_weighted_custom)
        
        gss_random_weighted_custom <- gss_random_weighted_custom %>%
        run_method(c(SS_SS_cssr # Stability selection (as proposed by Shah and
            # Samworth 2012)
            , SS_CSS_sparse_cssr # Sparse cluster stability selection
            # , SS_GSS_random_custom # Sparse cluster stability selection
            , SS_CSS_weighted_cssr # Weighted averaged cluster stability
            # selection
            , SS_CSS_avg_cssr # Simple averaged cluster stability
            # selection
            , clusRepLasso_cssr # Cluster representative lasso
            , protolasso_cssr # Protolasso
            , lasso_random # Lasso
            , elastic_net # Elastic Net
            , SS_CSS_sparse_cssr_est # Sparse cluster stability selection,
            # estimated clusters
            , SS_CSS_weighted_cssr_est # Weighted averaged cluster stability
            # selection, estimated clusters
            , SS_CSS_avg_cssr_est # Simple averaged cluster stability
            # selection, estimated clusters
            , clusRepLasso_cssr_est # Cluster representative lasso, estimated
            # clusters
            , protolasso_cssr_est # Protolasso, estimated clusters
        )) 

        save_simulation(gss_random_weighted_custom)

        gss_random_weighted_custom <- evaluate(gss_random_weighted_custom,
            list(cssr_mse))
    }


    save_simulation(gss_random_weighted_custom)
    

} else{
    print("loading simulation...")
    gss_random_weighted_custom <- load_simulation("gss_random_weighted_custom")
    print("simulation laoded!")
}

### Generate figures

results <- genPlotDf(gss_random_weighted_custom)

results_df <- results$results_df

n_methods <- results$n_methods

### Figure 4 (previously Figure 5) (known clusters)

fig_4_left <- createLossesPlot3(results_df[!(results_df$Method %in%
    nameMap(c("SS_CSS_sparse_cssr", "SS_CSS_avg_cssr", est_cluster_meths))), ],
    n_methods - 2 - length(est_cluster_meths))

fig_4_mid <- createNSBStabPlot2(results_df[!(results_df$Method %in%
    nameMap(c("SS_CSS_sparse_cssr", "SS_CSS_avg_cssr", est_cluster_meths))), ])

fig_4_right <- createStabMSEPlot2(results_df[!(results_df$Method %in%
    nameMap(c("SS_CSS_sparse_cssr", "SS_CSS_avg_cssr", est_cluster_meths))), ],
    n_methods - 2 - length(est_cluster_meths))

# 2. Save the legend
#+++++++++++++++++++++++
legend <- get_legend(fig_4_left + theme(legend.direction="horizontal"))

# 3. Remove the legend from the box plot
#+++++++++++++++++++++++
fig_4_left <- fig_4_left + theme(legend.position="none")

fig_4_mid <- fig_4_mid + theme(legend.position="none")

fig_4_right <- fig_4_right + theme(legend.position="none")

# 4. Arrange ggplot2 graphs with a specific width

fig_4 <- grid.arrange(fig_4_left, fig_4_mid, fig_4_right, legend, ncol=3,
    nrow = 2, layout_matrix = rbind(c(1, 2, 3), c(4, 4, 4)),
    widths = c(1.8, 1.8, 1.8), heights = c(2.5, 0.2))

fig_4 <- cowplot::ggdraw(fig_4) +
    theme(plot.background = element_rect(fill="white", color = NA))

print(fig_4)

saveFigure2(subdir="figures", plot=fig_4, size="large", filename="fig_4_known.pdf")

### Versions of Figure 4 plots with all methods (for supplement)

fig_4_supp_left <- createLossesPlot3(results_df[!(results_df$Method %in%
    nameMap(est_cluster_meths)), ], n_methods - length(est_cluster_meths))

saveFigure2(subdir="figures", plot=fig_4_supp_left, size="xmlarge",
    filename="sim_2_known_mse_supp.pdf")

fig_4_supp_mid <- createNSBStabPlot2(results_df[!(results_df$Method %in%
    nameMap(est_cluster_meths)), ])

saveFigure2(subdir="figures", plot=fig_4_supp_mid, size="xmlarge",
    filename="sim_2_known_stab_supp.pdf")

fig_4_supp_right <- createStabMSEPlot2(results_df[!(results_df$Method %in%
    nameMap(est_cluster_meths)), ], n_methods - length(est_cluster_meths))

saveFigure2(subdir="figures", plot=fig_4_supp_right, size="xmlarge",
    filename="sim_2_known_mse_stab_supp.pdf")



### Figure 4 (previously Figure 5) (estimated clusters)

fig_4_left <- createLossesPlot3(results_df[!(results_df$Method %in%
    nameMap(c("SS_CSS_sparse_cssr_est", "SS_CSS_avg_cssr_est",
        known_cluster_meths))), ], n_methods - 2 - length(known_cluster_meths))

fig_4_mid <- createNSBStabPlot2(results_df[!(results_df$Method %in%
    nameMap(c("SS_CSS_sparse_cssr_est", "SS_CSS_avg_cssr_est",
        known_cluster_meths))), ])

fig_4_right <- createStabMSEPlot2(results_df[!(results_df$Method %in%
    nameMap(c("SS_CSS_sparse_cssr_est", "SS_CSS_avg_cssr_est",
        known_cluster_meths))), ], n_methods - 2 - length(known_cluster_meths))

# 2. Save the legend
#+++++++++++++++++++++++
legend <- get_legend(fig_4_left + theme(legend.direction="horizontal"))

# 3. Remove the legend from the box plot
#+++++++++++++++++++++++
fig_4_left <- fig_4_left + theme(legend.position="none")

fig_4_mid <- fig_4_mid + theme(legend.position="none")

fig_4_right <- fig_4_right + theme(legend.position="none")

# 4. Arrange ggplot2 graphs with a specific width

fig_4 <- grid.arrange(fig_4_left, fig_4_mid, fig_4_right, legend, ncol=3,
    nrow = 2, layout_matrix = rbind(c(1, 2, 3), c(4, 4, 4)),
    widths = c(1.8, 1.8, 1.8), heights = c(2.5, 0.2))

fig_4 <- cowplot::ggdraw(fig_4) +
    theme(plot.background = element_rect(fill="white", color = NA))

print(fig_4)

saveFigure2(subdir="figures", plot=fig_4, size="large", filename="fig_4_est.pdf")

### Versions of Figure 4 plots with all methods (for supplement)

fig_4_supp_left <- createLossesPlot3(results_df[!(results_df$Method %in%
    nameMap(known_cluster_meths)), ], n_methods - length(known_cluster_meths))

saveFigure2(subdir="figures", plot=fig_4_supp_left, size="xmlarge",
    filename="sim_2_est_mse_supp.pdf")

fig_4_supp_mid <- createNSBStabPlot2(results_df[!(results_df$Method %in%
    nameMap(known_cluster_meths)), ])

saveFigure2(subdir="figures", plot=fig_4_supp_mid, size="xmlarge",
    filename="sim_2_est_stab_supp.pdf")

fig_4_supp_right <- createStabMSEPlot2(results_df[!(results_df$Method %in%
    nameMap(known_cluster_meths)), ], n_methods - length(known_cluster_meths))

saveFigure2(subdir="figures", plot=fig_4_supp_right, size="xmlarge",
    filename="sim_2_est_mse_stab_supp.pdf")

print("Total time:")

print(Sys.time() - t0)
