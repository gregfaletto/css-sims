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
n_sims <- 1000
# n_sims <- 100
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

    save_simulation(gss_random_ranking_custom_test0)
    

} else{
    gss_random_ranking_custom_test0 <- load_simulation("gss_random_ranking_custom_test0")
}

#### Generating figures

results <- genPlotDf(gss_random_ranking_custom_test0)

results_df <- results$results_df

n_methods <- results$n_methods

# Standalone proportion of subsamples figure

prop_fig <- createPhatPlot2(output(gss_random_ranking_custom_test0,
    methods="SS_SS_cssr")@out)

saveFigure2(subdir="figures", plot=prop_fig, size="slide",
    filename="ss_props.pdf")

### Figure 2

fig_2_right <- createLossesPlot3(results_df[results_df$Method %in%
    nameMap(c("SS_SS_cssr", "lasso_random")), ], 2)

# 2. Save the legend
#+++++++++++++++++++++++
legend_prop_fig <- get_legend(prop_fig +
    theme(legend.direction="horizontal"))

legend_fig_2_right <- get_legend(fig_2_right +
    theme(legend.direction="horizontal"))
# 3. Remove the legend from the box plot
#+++++++++++++++++++++++
fig_2_right_no_legend <- fig_2_right + theme(legend.position="none")
prop_fig_no_legend <- prop_fig + theme(legend.position="none")

# 4. Arrange ggplot2 graphs with a specific width

fig_2 <- grid.arrange(prop_fig_no_legend,
    fig_2_right_no_legend, legend_prop_fig, legend_fig_2_right
    , ncol=2, nrow=2, widths=c(0.5, 0.5), heights = c(2.5, 0.2)
    )

#  # 4. Arrange ggplot2 graphs with a specific width

fig_2 <- cowplot::ggdraw(fig_2) + theme(plot.background =
    element_rect(fill="white", color = NA))

saveFigure2(subdir="figures", plot=fig_2, size="mlarge", filename="fig_2.pdf")

### Figure 3

fig_3_left <- createLossesPlot3(results_df[!(results_df$Method %in%
    nameMap(c("SS_CSS_sparse_cssr", "SS_CSS_avg_cssr"))), ], n_methods - 2)

fig_3_mid <- createNSBStabPlot2(results_df[!(results_df$Method %in%
    nameMap(c("SS_CSS_sparse_cssr", "SS_CSS_avg_cssr"))), ])

fig_3_right <- createStabMSEPlot2(results_df[!(results_df$Method %in%
    nameMap(c("SS_CSS_sparse_cssr", "SS_CSS_avg_cssr"))), ], n_methods - 2)

# 2. Save the legend
#+++++++++++++++++++++++
legend <- get_legend(fig_3_left + theme(legend.direction="horizontal"))

# 3. Remove the legend from the box plot
#+++++++++++++++++++++++
fig_3_left <- fig_3_left + theme(legend.position="none")

fig_3_mid <- fig_3_mid + theme(legend.position="none")

fig_3_right <- fig_3_right + theme(legend.position="none")

# 4. Arrange ggplot2 graphs with a specific width

fig_3 <- grid.arrange(fig_3_left, fig_3_mid, fig_3_right, legend, ncol=3,
    nrow = 2, layout_matrix = rbind(c(1, 2, 3), c(4, 4, 4)),
    widths = c(1.8, 1.8, 1.8), heights = c(2.5, 0.2))

fig_3 <- cowplot::ggdraw(fig_3) +
    theme(plot.background = element_rect(fill="white", color = NA))

print(fig_3)

saveFigure2(subdir="figures", plot=fig_3, size="large", filename="fig_3.pdf")

### Versions of Figure 3 plots with all methods (for supplement)

fig_3_supp_left <- createLossesPlot3(results_df, n_methods)

saveFigure2(subdir="figures", plot=fig_3_supp_left, size="xmlarge",
    filename="sim_1_mse_supp.pdf")

fig_3_supp_mid <- createNSBStabPlot2(results_df)

saveFigure2(subdir="figures", plot=fig_3_supp_mid, size="xmlarge",
    filename="sim_1_stab_supp.pdf")

fig_3_supp_right <- createStabMSEPlot2(results_df, n_methods)

saveFigure2(subdir="figures", plot=fig_3_supp_right, size="xmlarge",
    filename="sim_1_mse_stab_supp.pdf")

# saveFigure()


# createLossesPlot3(results_df[results_df$Method %in% nameMap(c("SS_SS_cssr",
#     "lasso_random")), ], 2)

# createLossesPlot3(results_df, n_methods)

# createNSBStabPlot2(results_df)

# createStabMSEPlot2(results_df, n_methods)

# createPhatPlot2(output(completed_sim, methods="SS_SS_cssr")@out)



print("Total time:")

t1 <- Sys.time() - t0

print(t1)

print("Time per simulation:")

print(t1/n_sims)


