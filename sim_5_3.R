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
# library(stabs)
# library(lars)
library(glmnet)
library(digest)
library(knitr)
library(ggplot2)
library(doParallel)
library(ccaPP)
library(Metrics)
library(gridExtra)
library(irlba)
# library(Matrix)
library(doMC)
library(parallel)
library(dplyr)
library(cowplot)
# library(hash)
# library(digest)
# library(hashmap)
# library(CVXR) # For convex optimization in SC4 screening criterion evaluation;
## can eliminate if no longer using in that function


# PMA requires library "impute" which is no longer on CRAN. Instructions
# to download as listed on 
# http://www.bioconductor.org/packages/release/bioc/html/impute.html:

# ## try http:// if https:// URLs are not supported
# source("https://bioconductor.org/biocLite.R")
# biocLite("impute")

library(PMA)
library(e1071)
# library(pls)
# library(glmnetUtils)

############# Load Functions #################

wd <- getwd()
code_dir <- paste(wd, "Helper Functions", sep="/")
# setwd("/Users/gregfaletto/Google Drive/Data Science/R/USC/Stability Selection/toy_example/New method")
# source(file="stabsel_copy.R")
# sim_dir <- "/Users/gregfaletto/Google Drive/Data Science/R/USC/Stability Selection/toy_example"
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

# Folder in which to save figures
folder <- "210523/Weighted Sim (presentation slides)"

folder_main <- "/Users/gregfaletto/Dropbox/Subspace Stability Selection/sims_6"
folder_dir <- file.path(folder_main, folder)
dir.create(folder_dir, showWarnings = FALSE, recursive = TRUE)
# Run simulation, or load simulation that has been previously run?
run_new_sim <- FALSE
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
n_test <- 10000
n_sims <- 1000
# n_sims <- 2
# p <- 50
p <- as.list(c(0.5*n_model
    # , n_model
    # , 2*n_model 
    # , 5*n_model
    ))
nblocks <- 1
sig_blocks <- 1
k_unblocked <- 10
beta_low <- 1
# beta_high <- 2
beta_high <- as.list(c(
    # 1,
    1.5
    # , 2
    # , sqrt(k_unblocked)
    ))
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

# Display names of possible methods
DISPLAY_NAMES <- ("Lasso", "Stability Selection", "Sparse CSS",
    "Weighted Averaged CSS", "Simple Averaged CSS",
    "Cluster Representative Lasso", "Protolasso")

# if(!("lasso_random" %in% methods_to_eval)){
#     stop("!(lasso_random %in% methods_to_eval)")
# }

# Calculate clustered stability metric?

calc_clust_stab <- TRUE

if(calc_clust_stab & (sig_blocks != 1)){
    stop("sig_blocks != 1 (can't calculate clustered stability metric)")
}

setwd(sim_dir)

t0 <- Sys.time()

if(run_new_sim){

    print("Starting simulations...")

    set.seed(seed1)

    if(method=="MB"){
        gss_random_weighted_custom <- new_simulation("gss_random_weighted_custom",
        "GSS (Random Design, Weighted Averaging)") %>%
        generate_model(make_model=make_blocked_lin_mod4_ran_weight,
        n = n_model,
        p = p,
        # p = as.list(c(2*n_model, 40*n_model)),
        k_unblocked = k_unblocked,
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
        p = p,
        k_unblocked = k_unblocked,
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
        , vary_along = c("p", "beta_high")
        ) %>% simulate_from_model(nsim = n_sims) %>%
        run_method(c(lasso_random # lasso
            , SS_SS_cssr # Stability selection (as proposed by Shah and
            # Samworth 2012)
            , clusRepLasso_cssr
            # , BRVZ_avg_unwt
            , protolasso_cssr
            # , SS_CSS_sparse_cssr
            , SS_CSS_weighted_cssr
            # , SS_CSS_avg_cssr
            , elastic_net
        )) %>% evaluate(list(phat, labels))
    }


    save_simulation(gss_random_weighted_custom)
    

} else{
    gss_random_weighted_custom <- load_simulation("gss_random_weighted_custom")
}

completed_sim <- gss_random_weighted_custom

rm(gss_random_weighted_custom)

# Which simulation type is this? (Need for later processing)
sim <- "random_weighted"


# gss_random_weighted_custom <- run_method(gss_random_weighted_custom,
#     methods=c(lassoSS_phat_cor_squared)) %>% evaluate(list(phat, labels))




# # Getting n, p from original data generation
# n <- model(gss_random_weighted_custom)@params$n
# p <- model(gss_random_weighted_custom)@params$p

# gss_random_weighted_custom <- run_method(gss_random_weighted_custom, methods=c(lambda.max.cancor.ss)) 

# save_simulation(sim_6e)

# save_simulation(sim_6e)
# print("running evaluations...")
# sim_6e <- evaluate(sim_6e, list(sc4linear, size))
# save_simulation(sim_6e)


# if(n_param_combs == 1){
#     if(method=="MB"){
#         createPhatPlot(sim_model=model(gss_random_weighted_custom),
#             sim_eval=evals(gss_random_weighted_custom),
#             method="lassoMB_phat_random", ylab="Proportion of Subsamples",
#             n_iters=10, title=p_hat_titles)

#         createPhatPlot(sim_model=model(gss_random_weighted_custom),
#             sim_eval=evals(gss_random_weighted_custom),
#             method="lassoMB_phat_random", ylab="Proportion of Subsamples",
#             tau=0.51, line=TRUE, title=p_hat_titles)

#         # createPhatPlotLine(sim_model=model(gss_random_weighted_custom),
#         #     sim_eval=evals(gss_random_weighted_custom), method="lassoMB_phat",
#         #     ylab="Proportion of Subsamples", tau=0.51)

#         createPhatPlot(sim_model=model(gss_random_weighted_custom),
#             sim_eval=evals(gss_random_weighted_custom),
#             method="lassoMB_phat_ideal_random", ylab="Stability Score",
#             title=p_hat_titles)

#         # createPhatPlot(sim=gss_random_weighted_custom, method="lassoMB_phat_cor_squared",
#         #     ylab="Stability Score")
#     }else{
#         createPhatPlot(sim_model=model(gss_random_weighted_custom),
#             sim_eval=evals(gss_random_weighted_custom),
#             method="SS_SS_random", ylab="Proportion of Subsamples",
#             n_iters=10, title=p_hat_titles)

#         createPhatPlot(sim_model=model(gss_random_weighted_custom),
#             sim_eval=evals(gss_random_weighted_custom),
#             method="SS_SS_random", ylab="Proportion of Subsamples",
#             tau=0.51, line=TRUE, title=p_hat_titles)

#         # createPhatPlotLine(sim_model=model(gss_random_weighted_custom),
#         #     sim_eval=evals(gss_random_weighted_custom), method="lassoSS_phat",
#         #     ylab="Proportion of Subsamples", tau=0.51)

#         createPhatPlot(sim_model=model(gss_random_weighted_custom),
#             sim_eval=evals(gss_random_weighted_custom),
#             method="SS_GSS_random", ylab="Stability Score",
#             title=p_hat_titles)

#         # createPhatPlot(sim=gss_random_weighted_custom, method="lassoSS_phat_cor_squared",
#         #     ylab="Stability Score")
#     }
# }else{
#     for(p_prime in 1:n_param_combs){

#         # Get parameters

#         output_p_prime <- output(gss_random_weighted_custom)[[p_prime]]
#         params_p_prime <- model(gss_random_weighted_custom)[[p_prime]]@params
        
#         n_model_p_prime <- params_p_prime$n
#         p_p_prime <- params_p_prime$p
#         k_unblocked_p_prime <- params_p_prime$k_unblocked
#         beta_high_p_prime <- params_p_prime$beta_high
#         block_size_p_prime <- params_p_prime$block_size

#         # Skip this plot if doesen't match p_plot, beta_plot (if specified)

#         if(!is.na(p_plot)){
#             if(p_p_prime != p_plot){
#                 next
#             }
#         }

#         if(!is.na(beta_plot)){
#             if(beta_high_p_prime != beta_plot){
#                 next
#             }
#         }



#         if(method=="MB"){
#             createPhatPlot(sim_model=model(gss_random_weighted_custom)[[p_prime]],
#                 sim_eval=evals(gss_random_weighted_custom)[[p_prime]],
#                 method="lassoMB_phat_random", ylab="Proportion of Subsamples",
#                 n_iters=1, title=p_hat_titles)

#             # createPhatPlotLine(sim_model=model(gss_random_weighted_custom)[[p_prime]],
#             #     sim_eval=evals(gss_random_weighted_custom)[[p_prime]],
#             #     method="lassoMB_phat", ylab="Proportion of Subsamples", tau=0.51)

#             createPhatPlot(sim_model=model(gss_random_weighted_custom)[[p_prime]],
#                 sim_eval=evals(gss_random_weighted_custom)[[p_prime]],
#                 method="lassoMB_phat_random", ylab="Proportion of Subsamples",
#                 tau=0.51, line=TRUE, title=p_hat_titles)

#             createPhatPlot(sim_model=model(gss_random_weighted_custom)[[p_prime]],
#                 sim_eval=evals(gss_random_weighted_custom)[[p_prime]],
#                 method="lassoMB_phat_ideal_random", ylab="Stability Score",
#                 title=p_hat_titles)

#             # createPhatPlot(sim=gss_random_weighted_custom, method="lassoMB_phat_cor_squared",
#             #     ylab="Stability Score")
#         }else{
#             createPhatPlot(sim_model=model(gss_random_weighted_custom)[[p_prime]],
#                 sim_eval=evals(gss_random_weighted_custom)[[p_prime]],
#                 method="SS_SS_random", ylab="Proportion of Subsamples",
#                 n_iters=1, title=p_hat_titles)

#             # createPhatPlotLine(sim_model=model(gss_random_weighted_custom)[[p_prime]],
#             #     sim_eval=evals(gss_random_weighted_custom)[[p_prime]],
#             #     method="lassoSS_phat", ylab="Proportion of Subsamples", tau=0.51)

#             createPhatPlot(sim_model=model(gss_random_weighted_custom)[[p_prime]],
#                 sim_eval=evals(gss_random_weighted_custom)[[p_prime]],
#                 method="SS_SS_random", ylab="Proportion of Subsamples",
#                 tau=0.51, line=TRUE, title=p_hat_titles)

#             createPhatPlot(sim_model=model(gss_random_weighted_custom)[[p_prime]],
#                 sim_eval=evals(gss_random_weighted_custom)[[p_prime]],
#                 method="SS_GSS_random", ylab="Stability Score",
#                 title=p_hat_titles)

#             # createPhatPlot(sim=gss_random_weighted_custom, method="lassoSS_phat_cor_squared",
#             #     ylab="Stability Score")
#         }
#     }
# }



setwd(code_dir)

source("toy_example_plots.R")
