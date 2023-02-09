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
if(!is.null(dev.list())){
    dev.off()
}
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
library(gridExtra)
library(irlba)
library(Matrix)
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
registerDoMC(cores = detectCores() - 1)

setwd(code_dir)
source("toy_ex_slide_funcs.R")
source("model_functions.R")
source("method_functions.R")
source("eval_functions.R")
setwd(wd)



###############################################

# # Folder in which to save figures
# folder <- "210523/Ranking Sim (first simulation)"

# folder_main <- "/Users/gregfaletto/Dropbox/Subspace Stability Selection/sims_6"
# folder_dir <- file.path(folder_main, folder)
# dir.create(folder_dir, showWarnings = FALSE, recursive = TRUE)


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
n_test <- 1000
# n_sims <- 1000
n_sims <- 250
# p <- 50
p <- 100
# p <- as.list(c(0.5*n_model
#     # , n_model
#     # , 2*n_model 
#     # , 5*n_model
#     ))
nblocks <- 1
sig_blocks <- 1
k_unblocked <- 10
beta_low <- 1
# beta_low_min <- 0.55
# beta_low_max <- 1
# beta_high <- 2
beta_high <- 1.5
# beta_high <- as.list(c(
#     # 1,
#     1.5
#     # , 2
#     # , sqrt(k_unblocked)
#     ))
block_size <- 10
rho <- 0.9
# rho_low <- 0.9
# rho_high <- 0.9
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
        "SS_GSS_random_custom" # Sparse cluster stability selection
        , "SS_SS_random_custom" # Stability selection (as proposed by Shah and
        # Samworth 2012)
        , "lasso_random" # Lasso
        , "lasso_proto" # Protolasso
        , "SS_GSS_random_avg_unwt_custom" # Simple averaged cluster stability
        # selection
        , "BRVZ_avg_unwt" # Cluster representative lasso
        
        , "SS_GSS_random_avg_custom" # Weighted averaged cluster stability
        # selection
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
DISPLAY_NAMES <- c("Lasso", "Stability Selection", "Sparse CSS",
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

setwd(wd)

t0 <- Sys.time()

if(run_new_sim){

    print("Starting simulations...")

    set.seed(seed1)

    if(method=="MB"){
        gss_random_ranking_custom_test0 <- new_simulation("gss_random_ranking_custom_test0",
        "GSS (Random Design, ranked features)") %>%
        generate_model(make_model=make_sparse_blocked_linear_model4_random_ranking2,
        n = n_model,
        n_test = n_test,
        p = p,
        # p = as.list(c(2*n_model, 40*n_model)),
        k_unblocked = k_unblocked,
        beta_low = beta_low,
        beta_high = beta_high,
        nblocks = nblocks,
        sig_blocks = sig_blocks,
        block_size = block_size,
        rho = rho,
        # rho_low = rho_low,
        # rho_high = rho_high,
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
        ))

        gss_random_ranking_custom_test0 <- gss_random_ranking_custom_test0 %>%
            evaluate(list(cssr_mse))
            # evaluate(list(phat, labels))
    } else{
        gss_random_ranking_custom_test0 <- new_simulation("gss_random_ranking_custom_test0",
        "GSS (Random Design, ranked features)") %>%
        generate_model(make_model=make_sparse_blocked_linear_model4_random_ranking2,
        n = n_model,
        n_test = n_test,
        p = p,
        k_unblocked = k_unblocked,
        beta_low = beta_low,
        beta_high = beta_high,
        nblocks = nblocks,
        sig_blocks = sig_blocks,
        block_size = block_size,
        rho = rho,
        var = var,
        snr = snr 
        # , vary_along = ifelse(n_param_combs > 1, "p", NULL)
        # , vary_along = c("p", "beta_high")
        ) %>% simulate_from_model(nsim = n_sims) %>%
        run_method(c(
            SS_SS_cssr # Stability selection (as proposed by Shah and
            # Samworth 2012)
            # SS_SS_random_custom # Stability selection (as proposed by Shah and
            # Samworth 2012)
            , SS_CSS_sparse_cssr # Sparse cluster stability selection
            # , SS_GSS_random_custom # Sparse cluster stability selection
            , SS_CSS_weighted_cssr # Weighted averaged cluster stability
            # selection
            # , SS_GSS_random_avg_custom # Weighted averaged cluster stability
            # selection
            , SS_CSS_avg_cssr # Simple averaged cluster stability
            # selection
            # , SS_GSS_random_avg_unwt_custom # Simple averaged cluster stability
            # selection
            , clusRepLasso_cssr # Cluster representative lasso
            # , BRVZ_avg_unwt # Cluster representative lasso
            , protolasso_cssr # Protolasso
            # , lasso_proto # Protolasso
            , lasso_random # Lasso
            , elastic_net # Elastic Net
        ))

        gss_random_ranking_custom_test0 <- gss_random_ranking_custom_test0 %>%
            evaluate(list(cssr_mse))
            # evaluate(list(phat, labels))
    }

    # save_simulation(gss_random_ranking_custom_test0)
    

} else{
    gss_random_ranking_custom_test0 <- load_simulation("gss_random_ranking_custom_test0")
}

completed_sim <- gss_random_ranking_custom_test0

rm(gss_random_ranking_custom_test0)


# Which simulation type is this? (Need for later processing)
sim <- "random_ranking"

e <- evals(completed_sim)
edf <- as.data.frame(e)

methods <- unique(edf$Method)
n_methods <- length(methods)

alpha <- 0.05

model_sizes <- rep(1:(sig_blocks + k_unblocked), times=n_methods)
methods_vec <- nameMap(rep(methods, each=sig_blocks + k_unblocked))
mses <- rep(as.numeric(NA), (sig_blocks + k_unblocked)*n_methods)
margins <- rep(as.numeric(NA), (sig_blocks + k_unblocked)*n_methods)
nsbstabs <- rep(as.numeric(NA), (sig_blocks + k_unblocked)*n_methods)
nsb_lowers <- rep(as.numeric(NA), (sig_blocks + k_unblocked)*n_methods)
nsb_uppers <- rep(as.numeric(NA), (sig_blocks + k_unblocked)*n_methods)

# Get MSEs and NSB stabilities
for(i in 1:n_methods){
    edf_i <- edf[edf$Method == methods[i], ]
    # Should have a number of rows divisible by sig_blocks + k_unblocked
    stopifnot(nrow(edf_i) %% (sig_blocks + k_unblocked) == 0)
    stopifnot(n_sims*(sig_blocks + k_unblocked) == nrow(edf_i))

    meth_i_vec <- rep(as.numeric(NA), sig_blocks + k_unblocked)
    o_i <- output(completed_sim, methods=methods[i])@out

    for(k in 1:(sig_blocks + k_unblocked)){
        # MSE
        inds_k <- (sig_blocks + k_unblocked)*(0:(n_sims - 1)) + k
        stopifnot(all(inds_k %in% 1:nrow(edf_i)))
        mses_ik <- edf_i[inds_k, "MSE"]

        if(any(!is.na(mses_ik))){
            # mse_mat[k, i] <- mean(mses_ik, na.rm=TRUE)
            mses[(i - 1)*(sig_blocks + k_unblocked) + k] <- mean(mses_ik,
                na.rm=TRUE)
            # margin_mat[k, i] <- qnorm(1 - alpha/2)*sd(mses_ik, na.rm=TRUE)/
            #     sqrt(sum(!is.na(mses_ik)))
            margins[(i - 1)*(sig_blocks + k_unblocked) + k] <-
                qnorm(1 - alpha/2)*sd(mses_ik, na.rm=TRUE)/
                sqrt(sum(!is.na(mses_ik)))
        }

        # NSB Stability
        mat_i_k <- getBinMat(o_i, methods[i], k)
        # Get stability metric--only works if there are at least 2 simulations
        stopifnot(n_sims > 1)
        stab_res_ik <- calcNSBStabNone(mat_i_k, calc_errors=TRUE)
        nsbstabs[(i - 1)*(sig_blocks + k_unblocked) + k] <- stab_res_ik[1]
        nsb_lowers[(i - 1)*(sig_blocks + k_unblocked) + k] <- stab_res_ik[2]
        nsb_uppers[(i - 1)*(sig_blocks + k_unblocked) + k] <- stab_res_ik[3]
    }
}

results_df <- data.frame(ModelSize=model_sizes, Method=methods_vec, MSE=mses,
    MSELower=mses - margins, MSEUpper=mses + margins, NSBStability=nsbstabs,
    StabLower=nsb_lowers, StabUpper=nsb_uppers)

# createLossesPlot3(results_df, n_methods)

# createNSBStabPlot2(results_df)

# # NSB Stability
# results_list <- list()
# for(i in 1:n_methods){
#     # For each method and each model size, get a n_sims x p matrix of binary
#     # indicators of selected features (in matrix k, entry ij = 1 if feature j=
#     # was selected by the method on simulation i in the model of size k).
#     # This will be used as an input to the function calcNSBStab.
#     meth_i_vec <- rep(as.numeric(NA), sig_blocks + k_unblocked)
#     o_i <- output(completed_sim, methods=methods[i])@out
#     for(k in 1:(sig_blocks + k_unblocked)){
#         mat_i_k <- getBinMat(o_i, methods[i], k)
#         # Get stability metric--only works if there are at least 2 simulations
#         stopifnot(n_sims > 1)
#         meth_i_vec[k] <- calcNSBStabNone(mat_i_k)
#     }
#     results_list[[i]] <- meth_i_vec
# }

# names(results_list) <- methods

# saveFigure()

# createPhatPlot()

# createStabMSEPlot()


# # Getting n, p from original data generation
# n <- model(gss_random_ranking_custom_test0)@params$n
# p <- model(gss_random_ranking_custom_test0)@params$p

# setwd(code_dir)

# source("toy_example_plots.R")

print("Total time:")

t1 <- Sys.time() - t0

print(t1)

print("Time per simulation:")

print(t1/n_sims)


