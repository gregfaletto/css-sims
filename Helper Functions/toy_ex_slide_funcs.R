# setwd("/Users/gregfaletto/Google Drive/Data Science/LaTeX/Generalized Stability Selection Presentation")

library(ggplot2)

getClustersFromR <- function(R){
    # Check input

    stopifnot(is.matrix(R))
    stopifnot(all(dim(R) == p))
    stopifnot(all(diag(R) == 1))
    stopifnot(identical(R, t(R)))
    stopifnot(all(!is.na(R)))
    stopifnot(all(R %in% c(0, 1)))

    # Determine clusters from R
    clusters <- list()

    for(i in 1:nrow(R)){
        clusters[[i]] <- as.integer(which(R[i, ] > 0))
        stopifnot(length(clusters[[i]]) == length(unique(clusters[[i]])))
        stopifnot(all(!is.na(clusters[[i]])))
        stopifnot(is.integer(clusters[[i]]))
    }

    clusters <- unique(clusters)
    stopifnot(is.list(clusters))

    if(length(clusters) >= 2){
        # Check that clusters are non-overlapping
        for(i in 1:(length(clusters) - 1)){
            for(j in (i+1):length(clusters)){
                if(length(intersect(clusters[[i]], clusters[[j]])) != 0){
                    stop("Invalid R matrix with overlapping clusters (clusters must not be overlapping)")
                }
            }
        }
    }

    multiple <- FALSE

    stopifnot(is.list(clusters))
    # Only want to keep track of clusters with more than one member (we will
    # discard clusters with only one member)

    if(any(lengths(clusters) > 1)){ # & length(clusters) > 1
        
        # Only care about clusters with more than one element (only ones that need
        # to be treated differently)
        # keep track of whether there's more than one cluster or not
        multiple <- sum(lengths(clusters) > 1) > 1

        # Remove any clusters of size 1
        if(sum(lengths(clusters) > 1) == 1){
            clusters_temp <- clusters[[which(lengths(clusters) > 1)]]
            clusters <- list()
            clusters[[1]] <- clusters_temp
            rm(clusters_temp)
        } else if(multiple){
            clusters <- clusters[which(lengths(clusters) > 1)]
        }

        # Shouldn't be any clusters of size 1 at this point
        stopifnot(is.list(clusters))
        stopifnot(!multiple | all(lengths(clusters) > 1))
    } else{
        clusters <- list()
    }

    # Check output

    stopifnot(!multiple | all(lengths(clusters) > 1))
    stopifnot(is.list(clusters))
    stopifnot(all(lengths(clusters) > 1))
    stopifnot(all(!is.na(clusters)))
    stopifnot(length(clusters) == length(unique(clusters)))

    if(length(clusters) > 0){
        for(i in 1:length(clusters)){
            stopifnot(length(clusters[[i]]) == length(unique(clusters[[i]])))
            stopifnot(all(!is.na(clusters[[i]])))
            stopifnot(is.integer(clusters[[i]]) | is.numeric(clusters[[i]]))
            stopifnot(all(clusters[[i]] == round(clusters[[i]])))
        }

        if(length(clusters) >= 2){
            # Check that clusters are non-overlapping
            for(i in 1:(length(clusters) - 1)){
                for(j in (i+1):length(clusters)){
                    if(length(intersect(clusters[[i]], clusters[[j]])) != 0){
                        error_mes <- paste("Overlapping clusters detected; clusters must be non-overlapping. Overlapping clusters: ",
                            i, ", ", j, ".", sep="")
                        stop(error_mes)
                    }
                }
            }
        }
    } 

    return(clusters)
}

getTiedWithNext <- function(selected_clusts_names, clust_sel_props){
    p_selected <- length(selected_clusts_names)
    stopifnot(length(clust_sel_props) == p_selected)
    stopifnot(all(selected_clusts_names %in% names(clust_sel_props)))

    tied_with_next <- logical(p_selected)

    for(j in 1:(p_selected - 1)){
        final_sel_prop_j <- clust_sel_props[selected_clusts_names[j]]
        final_sel_prop_j_plus_1 <- clust_sel_props[selected_clusts_names[j + 1]]
        stopifnot(final_sel_prop_j >= final_sel_prop_j_plus_1)
        tied_with_next[j] <- final_sel_prop_j == final_sel_prop_j_plus_1
    }
    return(tied_with_next)
}

nameMap <- function(sys_name){
    # Takes in computer-friendly name for method and outputs display name

    # Inputs:

    # sys_name: A (vector of) character(s) with name(s) from method_functions.R

    # Outputs:

    # ret: A (vector of) character(s) with display friendly name(s)
    stopifnot(is.character(sys_name))
    ret <- sys_name
    ret[sys_name %in%  c("lasso", "lasso_random")] <- "Lasso"
    ret[sys_name %in%  c("lassoSS_phat", "SS_SS_random",
        "SS_SS_random_custom", "SS_SS", "SS_SS_cssr")] <- "Stability Selection"
    ret[sys_name %in%  c("lassoSS_phat_ideal", "SS_GSS_random",
        "SS_GSS_random_custom", "SS_GSS", "SS_CSS_sparse_cssr")] <- "Sparse CSS"
    ret[sys_name %in%  c("SS_CSS_sparse_cssr_est")] <- "Sparse CSS (est. clusts)"
    ret[sys_name %in%  c("SS_GSS_random_avg", "SS_GSS_random_avg_custom",
        "SS_GSS_avg", "SS_CSS_weighted_cssr")] <-
        "CSS"
        ret[sys_name %in%  c("SS_CSS_weighted_cssr_est")] <-
        "CSS (est. clusts)"
    ret[sys_name %in%  c("SS_GSS_random_avg_unwt", "SS_GSS_avg_unwt",
        "SS_GSS_random_avg_unwt_custom", "SS_CSS_avg_cssr")] <- "Simple Averaged CSS"
    ret[sys_name %in%  c("SS_CSS_avg_cssr_est")] <- "Simple Averaged CSS (est. clusts)"
    ret[sys_name %in%  c("BRVZ_avg_unwt", "clusRepLasso_cssr")] <-
        "Cluster Rep. Lasso"
        ret[sys_name %in%  c("clusRepLasso_cssr_est")] <-
        "Cluster Rep. Lasso (est. clusts)"
    ret[sys_name %in%  c("lasso_proto", "protolasso_cssr")] <- "Protolasso"
    ret[sys_name %in%  c("protolasso_cssr_est")] <- "Protolasso (est. clusts)"
    ret[sys_name %in%  c("elastic_net")] <- "Elastic Net"
    return(ret)
}

makeDfMethodsDisplay <- function(df){
    # Takes in data.frame with computer-friendly labels and returns 
    # data.frame with display-friendly labels

    stopifnot("Method" %in% colnames(df))

    labels <- as.factor(nameMap(as.character(df$Method)))
    labels <- droplevels(labels)
    df$Method <- labels

    stopifnot(all(levels(df$Method) %in% DISPLAY_NAMES))

    return(df)
}

createPhatPlot <- function(sim_model, sim_evals, method, ylab,
    feats_to_display=20, n_iters=1, line=FALSE, tau=NA, title=FALSE){
    # Creates bar plot of estimated probability of selection under our
    # proposed modification to stability selection against ranking of these
    # probabilities, in descending order of ranking (high to low estimated
    # probabilities). Bars are color coded according to whether feature
    # is blocked or unblocked, signal or noise.

    #
    # Inputs:
    #
    # sim_model: simulation model file with evaluations complete. The first simulation run
    # will be the one that is used for the plot; any further simulations will
    # be ignored.
    # sim_evals: simulation eval
    # method: character, method you would like to see plot for
    #
    # Outputs: ggplot2 object of plot

    require(ggplot2)

    if(line & is.na(tau)){
        stop("If line is TRUE, must provide tau")
    }

    p <- sim_model@params$p

    beta_high <- sim_model@params$beta_high

    if(n_iters > 1){
        # Create list of plots
        plots <- list()
    }

    for(i in 1:n_iters){
        # Index corresponding to this i
        id <- paste("r1.", i, sep="")

        data <- data.frame(1:p, sim_evals@evals[[method]][[id]]$phat,
            sim_evals@evals[[method]][[id]]$labels)

        colnames(data) <- c("Feature", "phat", "Legend")

        # re-name noise features

        levels(data$Legend) <- c(levels(data$Legend), "Noise Features",
            "Proxies For Z", "Weak Signal Features")

        data$Legend[data$Legend == "Unblocked Noise Features"] <- "Noise Features"
        data$Legend[data$Legend == "Blocked Noise Features"] <- "Proxies For Z"
        data$Legend[data$Legend == "Blocked Signal Features"] <- "Proxies For Z"
        data$Legend[data$Legend == "Unblocked Signal Features"] <- "Weak Signal Features"

        data <- drop(data)

        # Sort data in descending order of p_hat
        data <- data[order(data$phat, decreasing=TRUE), ]

        # Add a ranking variable
        colnames_orig <- colnames(data)
        data <- data.frame(1:p, data)
        colnames(data) <- c("Ranking", colnames_orig)

        # Keep only first feats_to_display features (easier to visualize)
        data <- data[1:min(nrow(data), feats_to_display), ]


        subtitle <- paste("p = ", p ,", beta_high = ", beta_high, sep="") 

        plot <- ggplot(data, aes(x=Ranking, y=phat, fill=Legend)) + geom_col() +
            ylab(ylab) 

        if(title){
            plot <- plot + labs(subtitle=subtitle)
        }

        if(line){
            plot <- plot + geom_hline(yintercept=tau, linetype="dashed",
                color="black") 
        }

        print(plot)

        if(n_iters > 1){
            # Create list of plots
            plots[[i]] <- plot
        }

    }
    if(n_iters > 1){
        return(plots)
    }
    return(plot)
}

createPhatPlotLine <- function(sim_model, sim_evals, method, ylab, tau,
    feats_to_display=20){
    # Creates bar plot of estimated probability of selection under our
    # proposed modification to stability selection against ranking of these
    # probabilities, in descending order of ranking (high to low estimated
    # probabilities). Bars are color coded according to whether feature
    # is blocked or unblocked, signal or noise.

    #
    # Inputs:
    #
    # sim_model: simulation model file with evaluations complete. The first simulation run
    # will be the one that is used for the plot; any further simulations will
    # be ignored.
    # sim_evals: simulation eval
    # method: character, method you would like to see plot for
    #
    # Outputs: ggplot2 object of plot

    require(ggplot2)

    p <- sim_model@params$p

    data <- data.frame(1:p, sim_evals@evals[[method]]$r1.1$phat,
        sim_evals@evals[[method]]$r1.1$labels)

    colnames(data) <- c("Feature", "phat", "Legend")

    # re-name noise features

    levels(data$Legend) <- c(levels(data$Legend), "Noise Features",
        "Proxies For Z", "Weak Signal Features")

    data$Legend[data$Legend == "Unblocked Noise Features"] <- "Noise Features"
    data$Legend[data$Legend == "Blocked Noise Features"] <- "Proxies For Z"
    data$Legend[data$Legend == "Blocked Signal Features"] <- "Proxies For Z"
    data$Legend[data$Legend == "Unblocked Signal Features"] <- "Weak Signal Features"

    data <- drop(data)

    # Sort data in descending order of p_hat
    data <- data[order(data$phat, decreasing=TRUE), ]

    # Add a ranking variable
    colnames_orig <- colnames(data)
    data <- data.frame(1:p, data)
    colnames(data) <- c("Ranking", colnames_orig)

    # Keep only first feats_to_display features (easier to visualize)
    data <- data[1:min(nrow(data), feats_to_display), ]
 
    plot <- ggplot(data, aes(x=Ranking, y=phat, fill=Legend)) + geom_col() +
        ylab(ylab) + geom_hline(yintercept=tau, linetype="dashed",
        color="black")
    print(plot)

    return(plot)
}

createPhatPlotTwoLines <- function(sim_model, sim_evals, method, ylab, tau1,
    tau2, feats_to_display=20){
    # Creates bar plot of estimated probability of selection under our
    # proposed modification to stability selection against ranking of these
    # probabilities, in descending order of ranking (high to low estimated
    # probabilities). Bars are color coded according to whether feature
    # is blocked or unblocked, signal or noise.
    
    #
    # Inputs:
    #
    # sim_model: simulation model file with evaluations complete. The first simulation run
    # will be the one that is used for the plot; any further simulations will
    # be ignored.
    # sim_evals: simulation eval
    # method: character, method you would like to see plot for
    #
    # Outputs: ggplot2 object of plot
    
    require(ggplot2)
    
    p <- sim_model@params$p
    
    data <- data.frame(1:p, sim_evals@evals[[method]]$r1.1$phat,
                       sim_evals@evals[[method]]$r1.1$labels)
    
    colnames(data) <- c("Feature", "phat", "Legend")
    
    # re-name noise features
    
    levels(data$Legend) <- c(levels(data$Legend), "Noise Features",
                            "Proxies For Z", "Weak Signal Features")
    
    data$Legend[data$Legend == "Unblocked Noise Features"] <- "Noise Features"
    data$Legend[data$Legend == "Blocked Noise Features"] <- "Proxies For Z"
    data$Legend[data$Legend == "Blocked Signal Features"] <- "Proxies For Z"
    data$Legend[data$Legend == "Unblocked Signal Features"] <- "Weak Signal Features"
    
    data <- drop(data)
    
    # Sort data in descending order of p_hat
    data <- data[order(data$phat, decreasing=TRUE), ]
    
    # Add a ranking variable
    colnames_orig <- colnames(data)
    data <- data.frame(1:p, data)
    colnames(data) <- c("Ranking", colnames_orig)

    # Keep only first feats_to_display features (easier to visualize)
    data <- data[1:min(nrow(data), feats_to_display), ]
    
    plot <- ggplot(data, aes(x=Ranking, y=phat, fill=Legend)) + geom_col() +
        ylab(ylab) + geom_hline(yintercept=tau1, linetype="dashed",
        color="black") + geom_hline(yintercept=tau2, linetype="dashed",
        color="black")
    print(plot)
    
    return(plot)
}

plotPredictionErrors <- function(sim, method){
    # Creates bar plot of estimated probability of selection under our
    # proposed modification to stability selection against ranking of these
    # probabilities, in descending order of ranking (high to low estimated
    # probabilities), along with prediction errors from model. 
    # Bars are color coded according to whether feature
    # is blocked or unblocked, signal or noise.

    #
    # Inputs:
    #
    # sim: simulation file with evaluations complete. The first simulation run
    # will be the one that is used for the plot; any further simulations will
    # be ignored.
    # method: character, method you would like to see plot for
    # 
    # Outputs: ggplot2 object of plot

    require(ggplot2)

    p <- model(sim)@params$p

    ev_6e <- evals(sim)

    data <- data.frame(1:p, ev_6e@evals[[method]]$r1.1$phat,
        ev_6e@evals[[method]]$r1.1$labels)

    colnames(data) <- c("Feature", "phat", "Legend")

    # Sort data in descending order of p_hat
    data <- data[order(data$phat, decreasing=TRUE), ]

    # Add a ranking variable
    colnames_orig <- colnames(data)
    data <- data.frame(1:p, data)
    colnames(data) <- c("Ranking", colnames_orig)

    # Add a column to data.frame for prediction errors
    colnames_orig <- colnames(data)
    data <- data.frame(data, numeric(p), numeric(p), numeric(p))
    colnames(data) <- c(colnames_orig, "PredictionError",
        "ScaledPredictionError", "ScreeningCriterion")

    # True model
    mu <- model(sim)@params$mu

    # Observed responses
    Y <- draws(sim)@draws$r1.1

    # For k in 1:p, calculate prediction error for top k variables as ranked by
    # stability selection relative to true model.
    for(k in 1:p){
        # Identify indices of top k features
        features_k <- data$Feature[1:k]
        X <- as.matrix(model(toy_model_new_meth_slide)@params$x[, features_k])
        # Predictions from OLS based on top k variables as ranked by stability
        # selection, using Y as response
        predictions_mu_k <- lm.fit(x=X, y=Y, tol=1e-12)$fitted.values
        # Calculate MSE (using mu as comparison)
        error_mu_k <- (1/n)*((predictions_mu_k - mu)%*%(predictions_mu_k - mu))
        data[k, "PredictionError"] <- error_mu_k
        data[k, "ScaledPredictionError"] <- error_mu_k/max(data$PredictionError)
        
        # Screening criterion of top k variables as ranked by stability
        # selection
        sc4_k <- sum(predictions_mu_k^2)/sum(mu^2)
        data[k, "ScreeningCriterion"] <- sc4_k
    }

    # Make labels display-friendly
    data <- makeDfLegendsDisplay(data)

    # Get max ranking value for integer labels on horizontal axis
    max_rank <- max(data$Ranking)
    if((max_rank %% 2) != 0){
        max_rank <- max_rank + 1
    }

    plot <- ggplot(data) + geom_col(mapping=aes(x=Ranking, y=phat, fill=Legend)) +
        geom_line(mapping=aes(x=Ranking, y=ScaledPredictionError)) +
        scale_x_continuous(breaks=seq(2, max_rank, by=2))

    print(plot)

    return(plot)
}



createFalseSelectionsPlot2 <- function(false_selecs, names, j_choices_mat,
    n_model, p, k, beta_high, J=NA, legend=TRUE){
    n_methods <- length(names)
    n_sims <- ncol(false_selecs)

    # If J is not provided, figure out on own
    stopifnot(length(J) <= 1)
    if(is.na(J)){
        J <- (nrow(false_selecs) - 1)/n_methods
    }
    stopifnot(J == round(J))
    stopifnot(all(false_selecs >= 0))
    stopifnot(!(n_sims != ncol(j_choices_mat) | nrow(false_selecs) != nrow(j_choices_mat)))

    if(n_sims > 1){
        row_inds <- which(rowSums(j_choices_mat) > 0)
        n_false_selecs <- sum(rowSums(j_choices_mat) > 0)
        false_selecs_vec <- numeric(n_false_selecs)
        j_choices <- which(rowSums(j_choices_mat[1:J, ]) > 0)
        stopifnot(length(j_choices)*n_methods == n_false_selecs)
        stopifnot(length(row_inds) == n_false_selecs)


        for(i in 1:n_false_selecs){
            cols_i <- j_choices_mat[row_inds[i], ]
            false_selecs_vec[i] <- mean(false_selecs[row_inds[i], cols_i])
        }
    } else{
        row_inds <- which(j_choices_mat)
        n_false_selecs <- sum(j_choices_mat)
        false_selecs_vec <- numeric(n_false_selecs)
        j_choices <- which(j_choices_mat[1:J, 1])
        stopifnot(length(j_choices)*n_methods == n_false_selecs)
        stopifnot(length(row_inds) == n_false_selecs)

        false_selecs_vec <- false_selecs[row_inds, 1]
    }

    df_gg <- data.frame(rep(j_choices, n_methods), false_selecs_vec,
        rep(names, each=length(j_choices)))

    # # Data.frame for plot
    # df_gg <- data.frame(rep(1:k_ss, n_ss_models), as.vector(false_selecs),
    #     rep(names, each=k))

    colnames(df_gg) <- c("Model_Size", "False_Selections", "Method")

    subtitle <- paste("n = ", n_model, ", p = ", p, ", k = ", k,
        ", beta_high = ", beta_high, sep="")

    # Make labels display-friendly
    df_gg <- makeDfMethodsDisplay(df_gg)

    # Get max ranking value for integer labels on horizontal axis
    max_rank <- max(df_gg$Model_Size)
    if((max_rank %% 2) != 0){
        max_rank <- max_rank + 1
    }

    plot <- ggplot(df_gg, aes(x=Model_Size, y=False_Selections, color=Method,
        shape=Method)) + suppressWarnings(geom_point(size=2.5, alpha=1)) + xlab("Model Size") +
        ylab("Number of False Selections") + labs(subtitle=subtitle) +
        scale_x_continuous(breaks=seq(2, max_rank, by=2))

    if(!legend){
        plot <- plot + theme(legend.position="none")
    }

    return(plot)
}



createProxySelectionsPlot2 <- function(proxy_selecs, names, j_choices_mat,
    n_model, p, k, beta_high, J=NA, legend=TRUE){
    n_methods <- length(names)
    n_sims <- ncol(proxy_selecs)

    # If J is not provided, figure out on own
    stopifnot(length(J) <= 1)
    if(is.na(J)){
        J <- (nrow(proxy_selecs) - 1)/n_methods
    }
    stopifnot(J == round(J))
    stopifnot(all(proxy_selecs >= 0))

    stopifnot(!(n_sims != ncol(j_choices_mat) | nrow(proxy_selecs) != nrow(j_choices_mat)))

    if(n_sims > 1){
        row_inds <- which(rowSums(j_choices_mat) > 0)
        n_proxy_selecs <- sum(rowSums(j_choices_mat) > 0)
        proxy_selecs_vec <- numeric(n_proxy_selecs)
        j_choices <- which(rowSums(j_choices_mat[1:J, ]) > 0)
        stopifnot(length(j_choices)*n_methods == n_proxy_selecs)

        stopifnot(length(row_inds) == n_proxy_selecs)


        for(i in 1:n_proxy_selecs){
            cols_i <- j_choices_mat[row_inds[i], ]
            proxy_selecs_vec[i] <- mean(proxy_selecs[row_inds[i], cols_i])
        }
    } else{
        row_inds <- which(j_choices_mat)
        n_proxy_selecs <- sum(j_choices_mat)
        proxy_selecs_vec <- numeric(n_proxy_selecs)
        j_choices <- which(j_choices_mat[1:J, 1])
        stopifnot(length(j_choices)*n_methods == n_proxy_selecs)
        stopifnot(length(row_inds) == n_proxy_selecs)

        proxy_selecs_vec <- proxy_selecs[row_inds, 1]
    }

    df_gg <- data.frame(rep(j_choices, n_methods), proxy_selecs_vec,
        rep(names, each=length(j_choices)))

    colnames(df_gg) <- c("Model_Size", "Proxy_Selections", "Method")

    # Make labels display-friendly
    df_gg <- makeDfMethodsDisplay(df_gg)

    # Get max ranking value for integer labels on horizontal axis
    max_rank <- max(df_gg$Model_Size)
    if((max_rank %% 2) != 0){
        max_rank <- max_rank + 1
    }

    subtitle <- paste("n = ", n_model, ", p = ", p, ", k = ", k, ", beta_high = ",
        beta_high, sep="")

    plot <- ggplot(df_gg, aes(x=Model_Size, y=Proxy_Selections, color=Method,
        shape=Method)) + suppressWarnings(geom_point(size=2.5, alpha=1)) + xlab("Model Size") +
        ylab("Number of Proxies Selected") + labs(subtitle=subtitle) +
        scale_x_continuous(breaks=seq(2, max_rank, by=2))

    if(!legend){
        plot <- plot + theme(legend.position="none")
    }

    return(plot)
}

createWeakSelectionsPlot2 <- function(proxy_selecs, false_selecs, names,
    j_choices_mat, n_model, p, k, beta_high, J=NA, legend=TRUE){
    n_methods <- length(names)
    n_sims <- ncol(proxy_selecs)

    # If J is not provided, figure out on own
    stopifnot(length(J) <= 1)
    if(is.na(J)){
        J <- (nrow(proxy_selecs) - 1)/n_methods
    }
    stopifnot(J == round(J))
    stopifnot(all(proxy_selecs >= 0))
    stopifnot(!(n_sims != ncol(j_choices_mat) | nrow(proxy_selecs) != nrow(j_choices_mat)))

    if(n_sims > 1){
        row_inds <- which(rowSums(j_choices_mat) > 0)
        n_proxy_selecs <- sum(rowSums(j_choices_mat) > 0)
        n_false_selecs <- sum(rowSums(j_choices_mat) > 0)
        proxy_selecs_vec <- numeric(n_proxy_selecs)
        false_selecs_vec <- numeric(n_false_selecs)
        j_choices <- which(rowSums(j_choices_mat[1:J, ]) > 0)
        stopifnot(length(j_choices)*n_methods == n_proxy_selecs)
        stopifnot(length(row_inds) == n_proxy_selecs)
        stopifnot(length(j_choices)*n_methods == n_false_selecs)
        stopifnot(length(row_inds) == n_false_selecs)

        for(i in 1:n_proxy_selecs){
            cols_i <- j_choices_mat[row_inds[i], ]
            proxy_selecs_vec[i] <- mean(proxy_selecs[row_inds[i], cols_i])
            false_selecs_vec[i] <- mean(false_selecs[row_inds[i], cols_i])
        }

    } else{
        row_inds <- which(j_choices_mat)
        n_proxy_selecs <- sum(j_choices_mat)
        n_false_selecs <- sum(j_choices_mat)
        proxy_selecs_vec <- numeric(n_proxy_selecs)
        false_selecs_vec <- numeric(n_false_selecs)
        j_choices <- which(j_choices_mat[1:J, 1])
        stopifnot(length(j_choices)*n_methods == n_proxy_selecs)
        stopifnot(length(row_inds) == n_proxy_selecs)
        stopifnot(length(j_choices)*n_methods == n_false_selecs)
        stopifnot(length(row_inds) == n_false_selecs)

        proxy_selecs_vec <- proxy_selecs[row_inds, 1]
        false_selecs_vec <- false_selecs[row_inds, 1]
        
    }

    # Now get number of weak signal selections: equal to model sizee minus
    # number of proxy selections minus number of false selections

    weak_selecs_vec <- rep(j_choices, n_methods) - proxy_selecs_vec -
        false_selecs_vec

    df_gg <- data.frame(rep(j_choices, n_methods), weak_selecs_vec,
        rep(names, each=length(j_choices)))

    # # Data.frame for plot
    # df_gg <- data.frame(rep(1:k_ss, n_ss_models), as.vector(proxy_selecs),
    #     rep(names, each=k))

    colnames(df_gg) <- c("Model_Size", "Weak_Signal_Selections", "Method")

    # Make labels display-friendly
    df_gg <- makeDfMethodsDisplay(df_gg)

    # Get max ranking value for integer labels on horizontal axis
    max_rank <- max(df_gg$Model_Size)
    if((max_rank %% 2) != 0){
        max_rank <- max_rank + 1
    }

    subtitle <- paste("n = ", n_model, ", p = ", p, ", k = ", k, ",
        beta_high = ", beta_high, sep="")

    plot <- ggplot(df_gg, aes(x=Model_Size, y=Weak_Signal_Selections,
        color=Method, shape=Method)) + suppressWarnings(geom_point(size=2.5, alpha=1)) +
        xlab("Model Size") + ylab("Number of Weak Signal Features Selected") + 
        labs(subtitle=subtitle) + scale_x_continuous(breaks=seq(2, max_rank,
            by=2))

    if(!legend){
        plot <- plot + theme(legend.position="none")
    }

    return(plot)

}

getErrorDf <- function(losses, names, j_choices_mat, J=NA,
    xaxis="No. Fitted Coefficients", cluster_count_mat=NA, plot_errors,
    conf=0.95){

    stopifnot(conf >= 0 & conf <= 1)

    n_methods <- length(names)
    n_sims <- ncol(losses)

    if(xaxis=="No. Fitted Coefficients"){

        names_vec <- rep(names, each=J)

        if(n_sims > 1){
            row_inds <- which(rowSums(j_choices_mat) > 0)
            n_losses <- sum(rowSums(j_choices_mat) > 0)
            losses_vec <- numeric(n_losses)
            if(plot_errors){
                lowers_vec <- numeric(n_losses)
                uppers_vec <- numeric(n_losses)
                alpha <- 1 - conf
            }
            j_choices <- which(rowSums(j_choices_mat) > 0)
            stopifnot(length(j_choices) == n_losses)
            stopifnot(length(row_inds) == n_losses)


            for(i in 1:n_losses){
                cols_i <- j_choices_mat[row_inds[i], ]
                losses_vec[i] <- mean(losses[row_inds[i], cols_i], na.rm=TRUE)
                if(plot_errors){
                    margin <- qnorm(1 - alpha/2)*sd(losses[row_inds[i], cols_i],
                        na.rm=TRUE)/sqrt(length(losses[row_inds[i], cols_i]))
                    lowers_vec[i] <- losses_vec[i] - margin
                    uppers_vec[i] <-  losses_vec[i] + margin
                }
            }
        } else{
            if(plot_errors){
                stop("can't have confidence intervals from loss if nsims=1")
            }
            row_inds <- which(j_choices_mat)
            n_losses <- sum(j_choices_mat)
            losses_vec <- numeric(n_losses)
            j_choices <- which(j_choices_mat[, 1])
            stopifnot(length(j_choices) == n_losses)
            stopifnot(length(row_inds) == n_losses)

            losses_vec <- losses[row_inds, 1]
            losses_vec <- losses_vec[!is.na(losses_vec)]
        }

        names_vec <- names_vec[row_inds]
        
        stopifnot(length(j_choices) == length(names_vec))
        stopifnot(length(j_choices) == length(losses_vec))

        no_coefs <- j_choices %% J
        no_coefs[no_coefs == 0] <- J

        df_gg <- data.frame(no_coefs, losses_vec, names_vec)

        if(plot_errors){
            df_gg <- data.frame(no_coefs, losses_vec, names_vec, lowers_vec,
                uppers_vec)
        }

    } else if(xaxis=="No. Selected Features"){
        if(plot_errors){
            stop("haven't written code to calculate confidence intervals for MSE when xaxis==No. Selected Features")
        }
        j_choices <- integer()
        losses_vec <- numeric()

        labels <- character()
        for(k in 1:n_methods){
            losses_vec_k <- as.vector(losses[((k - 1)*J + 1):(k*J), ])
            # Get all model sizes for this method
            model_sizes <- as.vector(cluster_count_mat[((k - 1)*J + 1):(k*J), ])

            stopifnot(length(losses_vec_k) == length(model_sizes))

            losses_vec_k <- losses_vec_k[model_sizes != 0]
            model_sizes <- model_sizes[model_sizes != 0]

            stopifnot(length(losses_vec_k) == length(model_sizes))

            model_sizes <- model_sizes[!is.na(losses_vec_k)]
            losses_vec_k <- losses_vec_k[!is.na(losses_vec_k)]
            
            stopifnot(length(losses_vec_k) == length(model_sizes))

            unique_model_sizes <- sort(unique(model_sizes))
   
            n_sizes_k <- length(unique_model_sizes)

            j_choices <- c(j_choices, unique_model_sizes)
            labels <- c(labels, rep(names[k], n_sizes_k))
            losses_means_k <- numeric(n_sizes_k)
            for(j in 1:n_sizes_k){
                size_j <- unique_model_sizes[j]
                stopifnot(sum(model_sizes == size_j) != 0)

                losses_means_k[j] <- mean(losses_vec_k[model_sizes == size_j],
                    na.rm=TRUE)

                stopifnot(!is.na(losses_means_k[j]))
            }
            losses_vec <- c(losses_vec, losses_means_k)
            stopifnot(!(length(j_choices) != length(losses_vec) | length(losses_vec) != length(labels)))
        }
        df_gg <- data.frame(j_choices, losses_vec, labels)

    }

    if(plot_errors){
        colnames(df_gg) <- c("Model_Size", "MSE", "Method", "MSELower", "MSEUpper")
    } else{
        colnames(df_gg) <- c("Model_Size", "MSE", "Method")
    }

    # Make labels display-friendly
    df_gg <- makeDfMethodsDisplay(df_gg)

    return(df_gg)
}


createLossesPlot2 <- function(losses, names, j_choices_mat, n_model=NA, p=NA,
    k=NA, beta_high=NA, J=NA, legend=TRUE, xaxis="No. Fitted Coefficients",
    cluster_count_mat=NA, plot_errors=FALSE, conf=0.95, subtitle=TRUE){

    n_methods <- length(names)
    n_sims <- ncol(losses)

    if(subtitle & any(is.na(c(n_model, p, k, beta_high)))){
        stop("if subtitle is TRUE, must provide n_model, p, k, beta_high")
    }

    stopifnot(conf >= 0 & conf <= 1)

    # If J is not provided, figure out on own
    stopifnot(length(J) <= 1)
    if(is.na(J)){
        J <- nrow(losses)/n_methods
    }
    stopifnot(J == round(J))
    stopifnot(!(any(losses[!is.na(losses)] < 0)))
    stopifnot(!(n_sims != ncol(j_choices_mat) | nrow(losses) != nrow(j_choices_mat)))

    if(xaxis=="No. Selected Features"){
        stopifnot(!(n_sims != ncol(cluster_count_mat) | nrow(losses) != nrow(cluster_count_mat)))
        stopifnot(all(cluster_count_mat[j_choices_mat] != 0))
    }

    if(xaxis=="No. Selected Features" & any(is.na(cluster_count_mat))){
        stop("If xaxis=No. Selected Features, must provide cluster_count_mat")
    }

    df_gg <- getErrorDf(losses, names, j_choices_mat, J, xaxis,
        cluster_count_mat, plot_errors=plot_errors, conf=conf)

    # Get max ranking value for integer labels on horizontal axis
    max_rank <- max(df_gg$Model_Size)
    if((max_rank %% 2) != 0){
        max_rank <- max_rank + 1
    }

    if(subtitle){
        subtitle_txt <- paste("n = ", n_model, ", p = ", p, ", k = ", k, ",
            beta_high = ", beta_high, sep="")
    }

    plot <- ggplot(df_gg, aes(x=Model_Size, y=MSE, color=Method,
        shape=Method)) + scale_shape_manual(values=1:n_methods) +
        suppressWarnings(geom_point(size=2.5, alpha=1)) + 
        # geom_hline(yintercept=true_model_loss
        #     # , color="red"
        #     , linetype = "dashed"
        #     ) + 
        xlab(xaxis) +
        scale_x_continuous(breaks=seq(2, max_rank, by=2))

    if(subtitle){
        plot <- plot + labs(subtitle=subtitle_txt)
    }

    if(plot_errors){
        plot <- plot + geom_errorbar(aes(ymin = MSELower, ymax = MSEUpper),
            width = 0.5)
    }

    if(!legend){
        plot <- plot + theme(legend.position="none")
    }

    return(plot)
}


createDesignForLoss <- function(X, y, selected, average=FALSE, to_avg=NA,
    avg_feats=NA, weights=NA){
    # Inputs:
    #
    # selected: a vector of length j containing a selected set.
    #
    # average: if TRUE, will do an averaging of cluster members.
    #

    # to_avg: must be provided if average=TRUE. A logical vector to_avg of the 
    # same length as 
    # selected. If feature j is in a cluster, the jth entry of to_avg will be
    # TRUE. 
    
    # avg_feats: must be provided if average=TRUE. A list of the same length as 
    # selected (or less) of integer vectors. 
    # If feature j is in a cluster, the jth entry of avg_feats will 
    # contain a vector of indices of the features from the
    # cluster of which feature j is a member.
    
    # weights: must be provided if average=TRUE. A list of the same length as 
    # avg_feats of numeric vectors.
    # If feature j is in a cluster, the jth entry of weights will
    # be the weights to use (in the same order as the jth entry of
    # avg_feats).

    p_selected <- length(selected)
    if(!average){
        X_df <- X[, selected]
    } else{
        X_df <- matrix(0, nrow(X), p_selected)
        colnames(X_df) <- character(p_selected)
        for(j in 1:p_selected){
            if(length(to_avg) >= j){
                if(to_avg[j]){
                    # Features to average
                    feats_to_avg_j <- avg_feats[[j]]
                    weights_j <- weights[[j]]
                    # Double-check that weights sum to 1
                    stopifnot(abs(sum(weights_j) - 1) < 10^(-6))
                    stopifnot(length(feats_to_avg_j) == length(weights_j))
                    if(length(weights_j) > 1){
                        X_df[, j] <- X[, feats_to_avg_j] %*% weights_j
                    } else{
                        stopifnot(length(weights_j) == 1)
                        X_df[, j] <- X[, feats_to_avg_j] * weights_j
                    }
                    stopifnot(length(avg_feats[[j]]) >= 1)
                    
                    if(p_selected > 1){
                        colnames(X_df)[j] <- paste("X", selected[j], sep="")
                    }
                } 
            } else{
                X_df[, j] <- X[, selected[j]]
            }
        }
    }

    df <- data.frame(y, X_df)
    colnames(df)[1] <- "y"

    # If only one feature was selected, then
    # X[, selected] is a vector, and the name of this column
    # of the data.frame will be weird. need to re-name the second column of
    # the data.frame in this case.
    if(p_selected == 1){
        colnames(df)[2] <- paste("X", selected, sep="")
    }

    return(df)
}

customMse <- function(actual, predicted){
    # Check inputs
    stopifnot(is.numeric(actual) | is.integer(actual))
    stopifnot(is.numeric(predicted) | is.integer(predicted))
    stopifnot(length(actual) == length(predicted))
    if(length(actual) == 0){
        return(0)
    }
    return(mean((actual - predicted)^2))
}

get_loss <- function(xtrain, ytrain, xtest, ytest, selected, average=FALSE,
    to_avg=NA, avg_feats=NA, weights=NA, digits=8) {
    # Inputs:
    #
    # selected: a vector of length j containing a selected set.
    #
    # average: if TRUE, will do an averaging of cluster members.
    #
    # to_avg: must be provided if average=TRUE. A logical vector to_avg of the 
    # same length as 
    # selected. If feature j is in a cluster, the jth entry of to_avg will be
    # TRUE. 
    
    # avg_feats: must be provided if average=TRUE. A list of the same length as 
    # selected (or less) of integer vectors. 
    # If feature j is in a cluster, the jth entry of avg_feats will 
    # contain a vector of indices of the features from the
    # cluster of which feature j is a member.
    
    # weights: must be provided if average=TRUE. A list of the same length as 
    # avg_feats of numeric vectors.
    # If feature j is in a cluster, the jth entry of weights will
    # be the weights to use (in the same order as the jth entry of
    # avg_feats).
    
    p_selected <- length(selected)
    if(p_selected == 0){
        return(NA)
    }
    n_train <- nrow(xtrain)
    n_test <- nrow(xtest)
    p <- ncol(xtrain)
    stopifnot(length(ytrain) == n_train)
    stopifnot(length(ytest) == n_test)
    stopifnot(p == ncol(xtest))
    stopifnot(max(selected) <= p)
    stopifnot(p_selected <= p)
    if(average){
        if(all(is.na(to_avg))){
            stop("If average==TRUE, must provide to_avg to get_loss")
        }
        if(any(to_avg)){
            if(all(is.na(avg_feats)) | all(is.na(weights))){
                stop("If average==TRUE, must provide avg_feats and weights to get_loss")
            }
        }
        stopifnot(length(to_avg) <= p)
        stopifnot(all.equal(lengths(avg_feats), lengths(weights)))
    }

    data_train <- createDesignForLoss(xtrain, ytrain, selected, average, to_avg,
        avg_feats, weights)

    data_test <- createDesignForLoss(xtest, ytest, selected, average, to_avg,
        avg_feats, weights)

    if(any(colnames(data_test)[2:ncol(data_test)] !=
        colnames(data_train)[2:ncol(data_train)])){
        colnames(data_test)[2:ncol(data_test)] <-
            colnames(data_train)[2:ncol(data_train)]
    }

    model <- lm(y~., data_train)
    loss <- customMse(actual=y_test, predicted=predict(model,
        newdata=data_test))

    return(round(loss, digits))
}

get_loss_mu <- function(xtrain, ytrain, mutrain, selected, average=FALSE,
    to_avg=NA, avg_feats=NA, weights=NA, digits=8){
    # Inputs:
    #
    # selected: a vector of length j containing a selected set.
    #
    # average: if TRUE, will do an averaging of cluster members.
    #
    # to_avg: must be provided if average=TRUE. A logical vector to_avg of the 
    # same length as 
    # selected. If feature j is in a cluster, the jth entry of to_avg will be
    # TRUE. 
    
    # avg_feats: must be provided if average=TRUE. A list of the same length 
    # (or less) as selected of integer vectors. 
    # If feature j is in a cluster, the jth entry of avg_feats will 
    # contain a vector of indices of the features from the
    # cluster of which feature j is a member.
    
    # weights: must be provided if average=TRUE. A list of the same length as 
    # avg_feats of numeric vectors.
    # If feature j is in a cluster, the jth entry of weights will
    # be the weights to use (in the same order as the jth entry of
    # avg_feats).
    
    p_selected <- length(selected)
    if(p_selected == 0){
        return(NA)
    }
    n_train <- nrow(xtrain)

    p <- ncol(xtrain)

    stopifnot(length(ytrain) == n_train)

    stopifnot(max(selected) <= p)
    stopifnot(p_selected <= p)
    if(average){
        if(length(to_avg) > 0 & all(is.na(to_avg))){
            stop("If average==TRUE, must provide to_avg to get_loss")
        }
        if(any(to_avg)){
            if((length(avg_feats) > 0 & all(is.na(avg_feats))) |
                (length(weights) > 0 & all(is.na(weights)))){
                stop("If average==TRUE, must provide avg_feats and weights to get_loss")
            }
        }
        stopifnot(length(to_avg) <= p)
        # print("lengths(avg_feats):")
        # print(lengths(avg_feats))
        # print("lengths(weights):")
        # print(lengths(weights))
        # for(i in 1:length(avg_feats)){
        #     if(lengths(avg_feats)[i] > 1){
        #         print("lengths(avg_feats)[i]:")
        #         print(lengths(avg_feats)[i])
        #         print("lengths(weights)[i]:")
        #         print(lengths(weights)[i])
        #         stopifnot(lengths(avg_feats)[i] == lengths(weights)[i])
        #     }
        # }
        # stopifnot(all.equal(unname(lengths(avg_feats)),
        #     unname(lengths(weights))))

        stopifnot(length(to_avg) <= p_selected)
        # stopifnot(length(weights) >= length(avg_feats))
        stopifnot(length(weights) <= p_selected)

        if(sum(to_avg) > 0){
            stopifnot(length(weights) >= max(which(to_avg)))
        }
    }

    data_train <- createDesignForLoss(xtrain, ytrain, selected, average, to_avg,
        avg_feats, weights)


    model <- lm(y~., data_train)
    loss <- customMse(actual=mutrain, predicted=predict(model))

    return(round(loss, digits))
}

num_false <- function(selected){
    if(length(selected) == 0){
        if(exists("wd")){
            setwd(wd)
        }
        stop("length(selected) == 0")
    }
    # Get correct selections
    num_cor_block_feats <- sig_blocks*block_size
    correct_selecs <- 1:(num_cor_block_feats + k_unblocked)
    # How many false selections?
    num_false <- sum(!(selected %in% correct_selecs))
    return(num_false)
}

num_proxies <- function(selected){
    if(length(selected) == 0){
        if(exists("wd")){
            setwd(wd)
        }
        stop("length(selected) == 0")
    }
    # Get proxy selections
    num_cor_block_feats <- sig_blocks*block_size
    correct_selecs <- 1:num_cor_block_feats
    # How many proxy selections?
    num_proxies <- sum(selected %in% correct_selecs)
    return(num_proxies)
}


identifyMethodIndex <- function(methodName, output){
    n_methods <- length(output)
    method_index <- 0
    for(i in 1:n_methods){
        if(output[[i]]@method_name==methodName){
            method_index <- i
        }
    }
    stopifnot(method_index != 0)
    return(method_index)
}


extractLassoSets <- function(selected, j_choices, j_choices_to_eliminate,
    var_names=NA){
    # List of selected sets to return
    selected_sets <- list()
    for(j in j_choices){
        # Lasso selected set of size j
        selected_j <- selected[lapply(selected, length) == j]
        if(length(selected_j) == 0){
            # First, check if we need to skip this j. True if the lasso path was unable to isolate a selected
            # set of size j (because the lambdas weren't dense enough--we could
            # re-run it in principle, but we'll just toss out that size set and
            # not worry about it.)
            # 
            # If we didn't find a lasso set of size j, just omit the stability
            # selected result and move on
            j_choices_to_eliminate <- c(j_choices_to_eliminate, j)
            next
        }
        selected_j <- selected_j[[1]]
        if(all(!is.na(var_names))){
            stopifnot(all(selected_j %in% 1:length(var_names)))
            stopifnot(max(selected_j) <= length(var_names))
            # Finally, record names of selected features (i.e., make this a  
            # vector of SNP names, not a vector of indices of features)
            selected_sets[[j]] <- var_names[selected_j]
        } else{
            selected_sets[[j]] <- selected_j
        }
        
    }
    return(list(selected_sets=selected_sets,
        j_choices_to_eliminate=j_choices_to_eliminate))
}

extractBRVZSets <- function(selected_sets, to_avg_list, selected_clusts_list,
    weights_list, j_choices, j_choices_to_eliminate, var_names=NA){

    selected_sets_ret <- list()
    to_avg_list_ret <- list()
    selected_clusts_list_ret <- list()
    weights_list_ret <- list()

    var_names_provided <- FALSE
    # If var_names is provided, convert avg_feats entries to character vectors
    if(all(!is.na(var_names))){
        stopifnot(all(!is.na(var_names)))
        var_names_provided <- TRUE
    }

    j_choices_to_eliminate <- c(j_choices_to_eliminate,
        j_choices[j_choices > max(vapply(selected_sets, length, integer(1)))])
    j_choices <- j_choices[j_choices <= max(vapply(selected_sets, length,
        integer(1)))]

    for(j in j_choices){
        if(length(selected_sets[[j]]) > 0){
            selected_sets_ret[[j]] <- selected_sets[[j]]

            selected_clusts_list_ret[[j]] <- selected_clusts_list[[j]]

            to_avg_list_ret[[j]] <- to_avg_list[[j]]
            
            weights_list_ret[[j]] <- weights_list[[j]]
        } else{
            j_choices_to_eliminate <- c(j_choices_to_eliminate, j)
        }
    }

    stopifnot(length(selected_sets) >= max(setdiff(j_choices,
        j_choices_to_eliminate)))

    return(list(selected_sets=selected_sets_ret,
        j_choices_to_eliminate=j_choices_to_eliminate,
        to_avg_list=to_avg_list_ret,
        selected_clusts_list=selected_clusts_list_ret, weights_list=weights_list_ret))
}

extractProtolassoSets <- function(selected_sets, j_choices,
    j_choices_to_eliminate, var_names=NA){

    selected_sets_ret <- list()

    var_names_provided <- FALSE
    # If var_names is provided, convert avg_feats entries to character vectors
    if(all(!is.na(var_names))){
        stopifnot(all(!is.na(var_names)))
        var_names_provided <- TRUE
    }

    j_choices_to_eliminate <- c(j_choices_to_eliminate,
        j_choices[j_choices > max(vapply(selected_sets, length, integer(1)))])
    j_choices <- j_choices[j_choices <= max(vapply(selected_sets, length,
        integer(1)))]

    for(j in j_choices){
        if(length(selected_sets[[j]]) > 0){
            selected_sets_ret[[j]] <- selected_sets[[j]]
        } else{
            j_choices_to_eliminate <- c(j_choices_to_eliminate, j)
        }
    }

    stopifnot(length(selected_sets_ret) <= max(setdiff(j_choices,
        j_choices_to_eliminate)))

    return(list(selected_sets=selected_sets_ret,
        j_choices_to_eliminate=j_choices_to_eliminate))
}

extractGSS_avgSets <- function(selected, to_avg, avg_feats, weights,
    tied_with_next, j_choices, j_choices_to_eliminate, var_names=NA){

    # Inputs

    # selected    
    # An integer vector of length at most p. Contains indices of features in 
    # decreasing order of weighted selection proportion. (Length may be less
    # than p because of clusters; only one member of each cluster will appear
    # in selected.)

    # to_avg: A logical vector to_avg of the 
    # same length as 
    # selected. If feature j is in a cluster, the jth entry of to_avg will be
    # TRUE. 

    # avg_feats: A list of the same length 
    # (or less) as selected of integer vectors. 
    # If feature j is in a cluster, the jth entry of avg_feats will 
    # contain a vector of indices of the features from the
    # cluster of which feature j is a member.
    
    # weights: A list of the same length as 
    # avg_feats of numeric vectors.
    # If feature j is in a cluster, the jth entry of weights will
    # be the weights to use (in the same order as the jth entry of
    # avg_feats).

    # tied_with_next
    # A logical vector the same length as selected. Entry j is TRUE if the jth
    # feature in selected has a tied selection proportion with feature j + 1
    # (i.e., no generalized stability selected model of size j exists).


    # Outputs:
    #
    # to_avg_list: A list whose elements are logical vectors of length j.
    # If feature k in 1:j is in a 
    # cluster, the kth entry of the jth vector in to_avg_list will be
    # TRUE. 
    #
    # selected_clusts_list: A list whose elements are lists of length j. If the model
    # of size j is defined, the jth entry of selected_clusts_list is a list containing
    # the first j list elements of avg_feats (with each entry converted to
    # character vectors containing the names of the corresponding variables
    # provided from var_names, if var_names is provided). That is, if feature k
    # in 1:j is in a cluster, the kth entry of the jth vector of selected_clusts_list
    # will contain (1) a character vector of the names of the features from the
    # cluster of which feature k is a member (if var_names is provided), or (2)
    # an integer vector of the indices of those features (if not).
    #
    # weights_list: A list whose elements are lists of length j. If
    # feature k in 
    # 1:j is in a cluster, the kth entry of the jth vector of
    # weights_list will
    # be the weights to use (in the same order as the jth entry of
    # selected_clusts_list).


    # Check inputs

    stopifnot(is.integer(selected))
    stopifnot(length(selected) == length(unique(selected)))

    stopifnot(is.logical(tied_with_next))
    stopifnot(length(tied_with_next) == length(selected))

    if(!is.null(to_avg)){
        stopifnot(is.logical(to_avg))
        stopifnot(length(to_avg) == length(selected))
        stopifnot(length(avg_feats) == length(weights))
        inds <- which(to_avg)
        stopifnot(is.list(avg_feats))
        stopifnot(is.list(weights))
    } else{
        inds <- integer()
    }

    for(i in inds){
        stopifnot(is.integer(avg_feats[[i]]))
        stopifnot(length(avg_feats[[i]]) >= 1)
        stopifnot(all(!is.na(avg_feats[[i]])))
        stopifnot(length(unique(avg_feats[[i]])) == length(avg_feats[[i]]))
        stopifnot((selected[i] %in% avg_feats[[i]]))
        stopifnot(is.numeric(weights[[i]]))
        stopifnot(length(weights[[i]]) >= 1)

        stopifnot(abs(sum(weights[[i]]) - 1) < 10^(-6))
    }

    var_names_provided <- FALSE
    # If var_names is provided, convert avg_feats entries to character vectors
    if(all(!is.na(var_names))){
        stopifnot(all(!is.na(var_names)))
        var_names_provided <- TRUE
        avg_feats_chars <- list()
        for(i in inds){
            avg_feats_chars[[i]] <- var_names[avg_feats[[i]]]
        }
        stopifnot(length(var_names) >= max(selected))
    }

    if(var_names_provided){
        stopifnot(is.list(avg_feats_chars))

        stopifnot(length(avg_feats_chars) == length(weights))

        for(i in inds){
            stopifnot(is.character(avg_feats_chars[[i]]))
            stopifnot(length(avg_feats_chars[[i]]) >= 2)
            stopifnot(all(!is.na(avg_feats_chars[[i]])))
            stopifnot(length(unique(avg_feats_chars[[i]])) == length(avg_feats_chars[[i]]))
            stopifnot(var_names[selected[i]] %in% avg_feats_chars[[i]])
        }
    }

    # List of selected sets to return
    selected_sets <- list()

    # Selected sets
    j_choices_to_eliminate <- c(j_choices_to_eliminate,
        j_choices[j_choices > length(selected)])
    j_choices <- j_choices[j_choices <= length(selected)]
    j_choices_to_eliminate <- c(j_choices_to_eliminate,
        which(tied_with_next))
    # Also need: to_avg, avg_feats, weights
    to_avg_list <- list()
    selected_clusts_list <- list()
    weights_list <- list()
    for(j in j_choices){
        selected_clusts_list[[j]] <- list()
        if(!tied_with_next[j]){
            if(var_names_provided){
                stopifnot(all(selected[1:j] %in% 1:length(var_names)))
                # Finally, record names of selected features (i.e., make this a  
                # vector of SNP names, not a vector of indices of features)
                selected_sets[[j]] <- var_names[selected[1:j]]


                # selected_clusts_list: A list whose elements are lists of length j. If
                # feature k in 
                # 1:j is in a cluster, the kth entry of the jth vector of selected_clusts_list
                # will contain the features from the
                # cluster of which feature k is a member.

                selected_clusts_list[[j]] <- avg_feats_chars[1:min(j,
                    length(avg_feats))]          

            } else{
                # The first j selected features
                selected_sets[[j]] <- selected[1:j]

                # selected_clusts_list: A list whose elements are lists of length j. If
                # feature k in 
                # 1:j is in a cluster, the kth entry of the jth vector of selected_clusts_list
                # will contain the features from the
                # cluster of which feature k is a member.

                selected_clusts_list[[j]] <- avg_feats[1:min(j, length(avg_feats))]
                
            }
            if(!is.null(to_avg)){
                # to_avg_list: A list whose elements are logical vectors of length j.
                # If feature k in 1:j is in a 
                # cluster, the kth entry of the jth vector in to_avg_list will be
                # TRUE. 
                to_avg_list[[j]] <- to_avg[1:min(j, length(to_avg))]
                
                # weights_list: A list whose elements are lists of length j. If
                # feature k in 
                # 1:j is in a cluster, the kth entry of the jth vector of
                # weights_list will
                # be the weights to use (in the same order as the jth entry of
                # selected_clusts_list).
                weights_list[[j]] <- weights[1:min(j, length(weights))]
            } else{
                to_avg_list[[j]] <- logical()
                weights_list[[j]] <- integer()
            }
        }
    }


    # Check output
    for(j in setdiff(j_choices, j_choices_to_eliminate)){

        if(!is.null(to_avg)){
            to_avg_j <- to_avg_list[[j]]
            avg_feats_j <- selected_clusts_list[[j]]
            weights_j <- weights_list[[j]]
            inds <- which(to_avg_j)
        } else{
            inds <- integer()
        }
        selected_j <- selected_sets[[j]]

        if(length(inds) > 0){
            stopifnot(length(avg_feats_j) >= max(inds))

            for(i in inds){
                if(var_names_provided){
                    stopifnot(is.character(avg_feats_j[[i]]))
                } else{
                    stopifnot(is.integer(avg_feats_j[[i]]))
                }
                
                stopifnot(length(avg_feats_j[[i]]) >= 1)

                stopifnot(all(!is.na(avg_feats_j[[i]])))

                stopifnot(length(unique(avg_feats_j[[i]])) == length(avg_feats_j[[i]]))
                stopifnot(selected_j[i] %in% avg_feats_j[[i]])
                stopifnot(is.numeric(weights_j[[i]]))
                stopifnot(length(weights_j[[i]]) >= 1)

                stopifnot(abs(sum(weights_j[[i]]) - 1) < 10^(-6))
            }
        }
    }

    return(list(selected_sets=selected_sets,
        j_choices_to_eliminate=j_choices_to_eliminate, to_avg_list=to_avg_list,
        selected_clusts_list=selected_clusts_list, weights_list=weights_list))
}

extractGSSSets <- function(selected, tied_with_next, j_choices,
    j_choices_to_eliminate, var_names=NA){

    selected_sets <- list()

    j_choices_to_eliminate <- c(j_choices_to_eliminate,
        j_choices[j_choices > length(selected)])
    j_choices <- j_choices[j_choices <= length(selected)]
    j_choices_to_eliminate <- c(j_choices_to_eliminate,
        which(tied_with_next))
    for(j in j_choices){
        if(!tied_with_next[j]){
            if(all(!is.na(var_names))){
                stopifnot(all(selected[1:j] %in% 1:length(var_names)))
                selected_sets[[j]] <- var_names[selected[1:j]]
            # The first j selected features
            } else{
                selected_sets[[j]] <- selected[1:j]
            }
        }
    }
    return(list(selected_sets=selected_sets,
        j_choices_to_eliminate=j_choices_to_eliminate))
}


extractSelectedSets <- function(out, method, j_choices){
    # Extracts selected sets given an @out object from simulator and
    # the name of the method

    # Outputs:

    # selected_sets
    # A list of vectors of selected sets. The jth entry of selected_sets
    # is an integer vector of length j, containing the indices of the j features
    # included in that selected set.

    # j_choices
    # A vector contaning the model sizes for which this method was defined
    # on this simulation run.
    
    # to_avg_list
    # A list whose elements are logical vectors of length j.
    # If feature k in 1:j is in a 
    # cluster, the kth entry of the jth vector in to_avg_list will be
    # TRUE. 
    #
    # selected_clusts_list
    # A list whose elements are lists of length j. If
    # feature k in 
    # 1:j is in a cluster, the kth entry of the jth vector of selected_clusts_list
    # will contain the features from the
    # cluster of which feature k is a member.
    #
    # weights_list
    # A list whose elements are lists of length j. If
    # feature k in 
    # 1:j is in a cluster, the kth entry of the jth vector of
    # weights_list will
    # be the weights to use (in the same order as the jth entry of
    # selected_clusts_list).

    # List of selected sets to return
    if(length(j_choices) == 0){
        warning("j_choices is empty")
        return(list(selected_sets=list(), j_choices=j_choices))
    }

    # J's to elminate from j_choices
    j_choices_to_eliminate <- integer()

    # Results from this iteration
    
    if(method %in% c("lasso", "lasso_random")){

        lasso_sets_results <- extractLassoSets(unique(out$selected), j_choices,
            j_choices_to_eliminate)

        selected_sets <- lasso_sets_results$selected_sets
        j_choices_to_eliminate <- lasso_sets_results$j_choices_to_eliminate

        rm(lasso_sets_results)

    } else if(method == "lasso_proto"){

        stopifnot(length(out$selected_sets) > 0)

        protolasso_results <- extractProtolassoSets(out$selected_sets, j_choices,
            j_choices_to_eliminate)

        selected_sets <- protolasso_results$selected_sets
        j_choices_to_eliminate <- protolasso_results$j_choices_to_eliminate

        rm(protolasso_results)

    } else if(method == "BRVZ_avg_unwt"){

        BRVZ_results <- extractBRVZSets(out$selected_sets, out$to_avg_list,
            out$selected_clusts_list, out$weights_list, j_choices,
            j_choices_to_eliminate)

        selected_sets <- BRVZ_results$selected_sets
        j_choices_to_eliminate <- BRVZ_results$j_choices_to_eliminate
        to_avg_list <- BRVZ_results$to_avg_list
        selected_clusts_list <- BRVZ_results$selected_clusts_list
        weights_list <- BRVZ_results$weights_list

        rm(BRVZ_results)

    } else if(method %in% c("SS_GSS_random_avg", # stability selection; no longer used
        "SS_GSS_random_avg_unwt", #no longer used
        "SS_GSS_random_avg_custom", # Weighted averaged cluster stability
        # selection
        "SS_GSS_random_avg_unwt_custom" # Simple averaged cluster stability
        )){

        GSS_avg_results <- extractGSS_avgSets(out$selected, out$to_avg,
            out$selected_clusts, out$weights, out$tied_with_next, j_choices,
            j_choices_to_eliminate)

        selected_sets <- GSS_avg_results$selected_sets
        j_choices_to_eliminate <- GSS_avg_results$j_choices_to_eliminate
        to_avg_list <- GSS_avg_results$to_avg_list
        selected_clusts_list <- GSS_avg_results$selected_clusts_list
        weights_list <- GSS_avg_results$weights_list

        rm(GSS_avg_results)
        
    } else if(method %in% c("SS_SS_random_custom", "SS_GSS_random_custom")){

        GSS_sets_results <- extractGSSSets(out$selected, out$tied_with_next,
            j_choices, j_choices_to_eliminate)

        selected_sets <- GSS_sets_results$selected_sets
        j_choices_to_eliminate <- GSS_sets_results$j_choices_to_eliminate

        rm(GSS_sets_results)
        
    } else if(method %in% c("SS_phat", "lassoMB_phat",
        "SS_SS_random")){

        ### Not using this currently

        selected_sets <- list()

        taus <- out$phat
        decreasing_taus <- sort(taus, decreasing=TRUE)
        # j_max <- sum(decreasing_taus > 0.5)
        # j_choices <- intersect(j_choices, 1:j_max)
        for(j in j_choices){

            # First, check if we need to skip this j. True if there is a tie in 
            # probabilities resulting in no stability selected set of size j
            # being defined.

            if(decreasing_taus[j] == decreasing_taus[j + 1]){
                # If there is a tie in probabilities, then there is no stability
                # selected set of size j defined. In that case, make a note of
                # this entry (will need to delete it later from both vectors) and
                #  move on to the next j.
                j_choices_to_eliminate <- c(j_choices_to_eliminate, j)
                next
            }

            # Stability selected set of size j    

            stability_selected_top_j_tau <- (decreasing_taus[j] +
                decreasing_taus[j + 1])/2

            stopifnot(stability_selected_top_j_tau < decreasing_taus[j] &
                stability_selected_top_j_tau > decreasing_taus[j + 1])

            selected_sets[[j]] <- which(taus > stability_selected_top_j_tau)

        }
    } else if(method %in% c("subspace_sel_grass", "subspace_ss_max_cancor",
        "subspace_ss_min_cancor")){

        ### Not using this currently
        selected_sets <- list()


        selected_set <- out$selected
        j <- length(selected_set)
        j_choices <- intersect(j_choices, j)
        if(length(j_choices) == 0){
            warning("j_choices is empty")
            return(list(selected_sets, j_choices))
        }
        selected_sets[[j]] <- selected_set
        # print(paste("q for ", method, ": ", out$q, sep=""))
    } else if(method %in% c("lassoSS_phat_ideal", "SS_GSS_random",
        "lassoMB_phat_ideal")){


        ### Not using this currently

        selected_sets <- list()

        # taus <- out$phat
        # if(sig_blocks != 1){
        #     stop("sig_blocks != 1")
        # }
        # if(nblocks > 1){
        #     stop("nblocks > 1")
        # }
        # p_feat <- length(taus)
        # if(nrow(out$R) != p_feat | ncol(out$R) != p_feat){
        #     stop("nrow(out$R) != p_feat | ncol(out$R) != p_feat")
        # }
        # names(taus) <- 1:p_feat
        # # Recover block structure: find features in a cluster of size
        # # more than 1. Assume only one cluster in data.
        # clustered_features <- numeric()
        # for(j in 1:(p_feat - 1)){
        #     for(k in (j+1):p_feat){
        #         if(out$R[j, k] != 0){
        #             clustered_features <- c(clustered_features, j, k)
        #         }
        #     }
        # }
        # clustered_features <- unique(clustered_features)
        # cluster_rep <- clustered_features[1]
        # # Need to create new kind of decreasing taus. Any features from
        # # cluster will count as one feature, and the tau will be capped at 1.
        # taus <- sort(taus, decreasing=TRUE)
        # features_ordered <- as.integer(names(taus))
        # # Remove clustered features (other than cluster_rep)
        # features_ordered <- setdiff(features_ordered,
        #     clustered_features[2: length(clustered_features)])
        # # Selected sets
        # for(j in j_choices){
        #     selected_sets[[j]] <- sort(features_ordered[1:j],
        #         decreasing=FALSE)
        # }

        # Modified code: selected from function output now has
        # features in decreasing order of weighted selection proportion
        # by default
        # Selected sets
        selected <- out$selected
        j_choices <- j_choices[j_choices <= length(selected)]
        for(j in j_choices){
            selected_sets[[j]] <- selected[1:j]
        }
        # print(paste("selected_sets for", nameMap(method), ":"))
        # print(selected_sets)
    } else{
        message <- paste("No feature extraction method exists for method",
            method, ".")
        stop(message)
    }


    if(length(j_choices_to_eliminate) > 0){
        j_choices <- setdiff(j_choices, j_choices_to_eliminate)
    } 

    if(all(is.null(selected_sets))){
        print("method:")
        print(method)
        stop("all(is.null(selected_results[[1]]))")
    }

    stopifnot(length(selected_sets) > 0)
    stopifnot(length(j_choices) > 0)

    if(!(method %in% c("SS_GSS_random_avg", "SS_GSS_random_avg_unwt",
        "BRVZ_avg_unwt", "SS_GSS_random_avg_custom",
        "SS_GSS_random_avg_unwt_custom"))){
        return(list(selected_sets=selected_sets, j_choices=j_choices))
    }

    return(list(selected_sets=selected_sets, j_choices=j_choices,
        to_avg_list=to_avg_list, selected_clusts_list=selected_clusts_list,
        weights_list=weights_list))
}


calcHamStab <- function(sel_mat, n_sims=NA){
    # Calculates the (empirical) instability of a feature selection method
    # from a matrix containing several selected sets. Does this by calculating
    # the Hamming distance between every pair of selected sets (that is, the
    # number of features where there is disagreement between the two selected
    # sets). If this metric equals p, then none of the pairs of selected sets
    # had any features in common; if this metric equals 0, every selected set
    # was identical. 

    # Inputs:

    # sel_mat: A matrix of sets of selected features. The number of columns
    # is equal to the number of features. Each row contains a selected set.
    # In selected set i, if feature j was selected then sel_mat[i, j] = 1 and
    # if feature j was not selected then sel_mat[i, j] = 0.

    # n_sims: Optional; should match number of rows of sel_mat. (Can provide
    # if desired for safety to confirm that nrow(sel_mat) == n_sims.) 

    # Returns:

    # The (empirical) Hamming instability of the feature selection method.
    # Ranges between 0 (complete stability of method; perfect agreement between
    # all selected sets) and p (all selected sets differed from each other).
    # Returns NA if unable to calculate metric (for example, if all selected
    # sets were empty).

    require(e1071)

    # Check inputs

    stopifnot(is.matrix(sel_mat))
    stopifnot(all(!is.na(sel_mat)))
    stopifnot(all(sel_mat %in% c(0, 1)))

    if(is.na(n_sims)){
        n_sims <- nrow(sel_mat)
    }

    stopifnot(nrow(sel_mat) == n_sims)
    p <- ncol(sel_mat)
    if(all(sel_mat == 0)){
        return(NA)
    }
    # Eliminate rows of all zeroes (no selections)
    all_zero_rows <- apply(sel_mat, 1, function(x){all(x == 0)})
    if(sum(!all_zero_rows) == 1){
        return(NA)
    }
    sel_mat_nonzero <- sel_mat[!all_zero_rows, ]
    stopifnot(ncol(sel_mat_nonzero) == p)
    m <- nrow(sel_mat_nonzero)
    stopifnot(m > 1)
    if(m <= n_sims/2){
        print("Warning here")
        warning("m <= n_sims/2")
    }
    ham_dist_mat <- e1071::hamming.distance(sel_mat_nonzero)

    ham_instab <- sum(ham_dist_mat[lower.tri(ham_dist_mat)])/(m*(m-1)/2)

    # Check output
    stopifnot(length(ham_instab) == 1)
    stopifnot(is.numeric(ham_instab) | is.integer(ham_instab))
    stopifnot(ham_instab >= 0)
    stopifnot(ham_instab <= p)

    return(ham_instab)
}

calcNSBStab <- function(sel_mat, n_sims=NA, calc_errors=FALSE, conf=0.95,
    coarseness=1, cluster_count="none", n_clust_feats=0, C=NA){
    # Calculates the (empirical) stability of a feature selection method
    # from a matrix containing several selected sets, using either the metric from
    # Nogueira et. al (2018) or, if accounting for clustered features is
    # desired, the method from Sechidis et. al (2019). In this metric, a 
    # completely stable feature selection algorithm has stability 1 (in 
    # practice, this happens when
    # every set selected by the algorithm from separate draws of the data is
    # identical). A feature selection algorithm that chooses features completely
    # at random has expected value 0 under this metric. (Negative values of the
    # metric are possible.) Optionally also calculates a confidence interval
    # for the stability metric.

    # References: 

    # Nogueira, S., Sechidis, K., & Brown, G. (2018). On the stability of
    # feature selection algorithms. Journal of Machine Learning Research, 18, 
    # #1–54. Retrieved from http://jmlr.org/papers/v18/17-514.html.

    # Sechidis, K., Papangelou, K., Nogueira, S., Weatherall, J. and Brown, G.,
    # 2019, June. On the Stability of Feature Selection in the Presence of
    # Feature Correlations. In ECML/PKDD (1) (pp. 327-342).

    # Inputs:

    # sel_mat: A matrix of sets of selected features. The number of columns
    # is equal to the number of features. Each row contains a selected set.
    # In selected set i, if feature j was selected then sel_mat[i, j] = 1 and
    # if feature j was not selected then sel_mat[i, j] = 0.

    # n_sims: TODO: check if this argument is no longer needed?

    # calc_errors: logical; whether or not to calculate a confidence interval
    # for the stability metric. Only possible if cluster_count is not "corr".

    # conf: numeric; if calc_errors is TRUE, then a confidence interval will
    # be calculated at this confidence level. (Must be strictly between 0 and
    # 1.)

    # coarseness: TODO: check if this argument is no longer needed?

    # cluster_count: character; one of "none", "single", "multi", or "corr".
    # Whether to account for clustering in the features when calculating the
    # stability metric (using the metric from Sechidis et. al 2019) or not
    # (using the metric from Nogueira et. al. 2018). If "none", disregards
    # clustering in features and uses metric from Nogueira et. al (2018). If
    # "corr", calculates stability metric from Sechidis et. al (2019). This
    # requires providing C (see below). "single" and "multi" are experimental
    # settings not justified by theory. For "single" and "multi", the metric 
    # assumes that the data contains one cluster of features represented by the
    # first n_clust_feats columns of sel_mat. If cluster_count is "single", then
    # for each row  of sel_mat if any of the first n_clust_feats features are 
    # selected,
    # the first feature in the cluster will be marked as selected (labeled 1)
    # and the rest will be labeled as 0. If cluster_count is "multi", then
    # for each row if k of the the first n_clust_feats features are selected
    # then features 1:k of sel_mat will be labeled as 1 and the rest will be
    # labeled as 0.

    # n_clust_feats: Required if and only if cluster_count is "single" or
    # "multi". Then the first n_clust_feats columns of sel_mat are assumed
    # to represent features in the same cluster. (See description of
    # cluster_count.)

    # C: Required if and only if cluster_count=="corr". A p x p matrix (where p
    # = ncol(sel_mat) is the  number of features) with all nonnegative entries 
    # between 0 and 1
    # containing the similarity of features. The diagonal entries should all
    # equal 1, and for i != j, C[i, j] is the similarity of feature i and j
    # (for example, one could choose C[i, j] to be the absolute value of the
    # correlation between features i and j).

    # Returns:

    # If calc_errors is FALSE (or cluster_count is "corr"): returns the value
    # of the stability metric, which is at most 1. (Returns NA if metric
    # is unable to be calculated, for example if sel_mat contains all 0s.) If
    # calc_errors is TRUE (and cluster_count is not "corr"), returns a length
    # 3 numeric vector containing the calculated stability metric, a lower
    # bound for the conf confidence interval, and an upper bound for the
    # confidence interval.

    if(is.na(n_sims)){
        n_sims <- nrow(sel_mat)/coarseness
    }

    # Check inputs
    checkNSBStabInputs(sel_mat, n_sims, calc_errors, conf, coarseness,
        cluster_count, n_clust_feats, C)
    p <- ncol(sel_mat)
    
    if(all(sel_mat == 0)){
        return(NA)
    }
    # Eliminate rows of all zeroes (no selections)
    all_zero_rows <- apply(sel_mat, 1, function(x){all(x == 0)})
    M <- sum(!all_zero_rows)
    if(M <= 1){
        return(NA)
    }

    if(M <= n_sims*coarseness/10){
        print("Warning here")
        warning("sum(!all_zero_rows) <= n_sims*coarseness/10")
    }
    sel_mat_nonzero <- sel_mat[!all_zero_rows, ]
    stopifnot(ncol(sel_mat_nonzero) == p)
    stopifnot(M == nrow(sel_mat_nonzero))

    if(cluster_count != "corr"){
        if(cluster_count %in% c("multi", "single")){
            # If modifying criterion for clusters: for each row, count how
            # many clustered features were selected
            n_clust_selecs <- rowSums(sel_mat_nonzero[, 1:n_clust_feats])
            sel_mat_nonzero[, 1:n_clust_feats] <- 0
            if(cluster_count == "multi"){
                for(i in 1:nrow(sel_mat_nonzero)){
                    # For the ith row, fill the first n_clust_selecs[i] cluster
                    # features with ones (i.e., don't distinguish amongst clustered
                    # features)
                    sel_mat_nonzero[i, 1:n_clust_selecs[i]] <- 1
                }
            } else if(cluster_count=="single"){
                for(i in 1:nrow(sel_mat_nonzero)){
                    # For the ith row, add a 1 to only the first cluster feature 
                    # (i.e., don't distinguish amongst clustered features)
                    sel_mat_nonzero[i, 1] <- 1
                }
                sel_mat_nonzero <- cbind(sel_mat_nonzero[, 1],
                    sel_mat_nonzero[,
                    (n_clust_feats + 1):ncol(sel_mat_nonzero)])
                p <- ncol(sel_mat)
            } else{
                stop("cluster_count must be single, multi, or none")
            }
        }

        # Calculate s_f_squared values
        p_hat <- colMeans(sel_mat_nonzero)
        s_f_squared <- apply(sel_mat_nonzero, 2,
            function(x){M/(M-1)*mean(x)*(1-mean(x))})
        stopifnot(all.equal(s_f_squared, M/(M-1)*p_hat*(1-p_hat)))
        k <- rowSums(sel_mat_nonzero)
        k_bar <- mean(k)

        stat <- 1 - mean(s_f_squared)/(k_bar/p*(1 - k_bar/p))

        # Check output
        stopifnot(length(stat) == 1)
        stopifnot(is.numeric(stat) | is.integer(stat))
        stopifnot(stat <= 1)

        if(calc_errors){
            alpha <- 1 - conf
            stopifnot(alpha >= 0 & alpha <= 1)

            # Calculate phi_hat_i's
            phi_hat_i <- numeric(M)
            stopifnot(ncol(sel_mat_nonzero) == length(p_hat))
            for(i in 1:M){
                phi_hat_i[i] <- 1/(k_bar/p*(1 - k_bar/p))*
                    (mean(sel_mat_nonzero[i, ]*p_hat) - k[i]*k_bar/p^2 +
                    stat/2*(2*k_bar*k[i]/p^2 - k[i]/p - k_bar/p + 1))
            }
            stopifnot(all(phi_hat_i != 0))
            v <- 4/M^2*sum((phi_hat_i - mean(phi_hat_i))^2)
            stopifnot(v >= 0)
            # Calculate errors for error bars
            margin <- qnorm(1 - alpha/2)*sqrt(v)

            # Check output

            stopifnot(length(margin) == 1)
            stopifnot(is.numeric(margin) | is.integer(margin))
            stopifnot(margin <= 1)
            stopifnot(margin >= 0)

            # Upper margin can't be greater than 1 since statistic can't be
            # greater than 1
            return(c(stat, stat - margin, min(stat + margin, 1)))
        }

        return(stat)
    } else if(cluster_count=="corr"){
        # Get k_bar
        k <- rowSums(sel_mat_nonzero)
        k_bar <- mean(k)
        # Calculate sample covariance
        S <- cov(sel_mat_nonzero)
        # Calculate Sigma_0

        # Off-diagonal entries
        off_diags <- (k_bar^2 - k_bar)/(p^2 - p) - k_bar^2/p^2
        # Diagonal entries
        diags <- k_bar/p*(1 - k_bar/p)

        Sigma_0 <- matrix(off_diags, p, p)
        diag(Sigma_0) <- diags

        # Stability metric
        stat <- 1 - sum(diag(C %*% S))/sum(diag(C %*% Sigma_0))
        
        # Check output
        stopifnot(length(stat) == 1)
        stopifnot(is.numeric(stat) | is.integer(stat))
        stopifnot(stat <= 1)

        return(stat)
    } else{
        stop("!(cluster_count %in% c(none, single, multi, corr))")
    }
    stop("Function calcNSBStab did not return value")
}

calcNSBStabNone <- function(sel_mat, n_sims=NA, calc_errors=FALSE, conf=0.95,
    coarseness=1){
    # Calculates the (empirical) stability of a feature selection method
    # from a matrix containing several selected sets, using either the metric from
    # Nogueira et. al (2018) or, if accounting for clustered features is
    # desired, the method from Sechidis et. al (2019). In this metric, a 
    # completely stable feature selection algorithm has stability 1 (in 
    # practice, this happens when
    # every set selected by the algorithm from separate draws of the data is
    # identical). A feature selection algorithm that chooses features completely
    # at random has expected value 0 under this metric. (Negative values of the
    # metric are possible.) Optionally also calculates a confidence interval
    # for the stability metric.

    # References: 

    # Nogueira, S., Sechidis, K., & Brown, G. (2018). On the stability of
    # feature selection algorithms. Journal of Machine Learning Research, 18, 
    # #1–54. Retrieved from http://jmlr.org/papers/v18/17-514.html.

    # Sechidis, K., Papangelou, K., Nogueira, S., Weatherall, J. and Brown, G.,
    # 2019, June. On the Stability of Feature Selection in the Presence of
    # Feature Correlations. In ECML/PKDD (1) (pp. 327-342).

    # Inputs:

    # sel_mat: A matrix of sets of selected features. The number of columns
    # is equal to the number of features. Each row contains a selected set.
    # In selected set i, if feature j was selected then sel_mat[i, j] = 1 and
    # if feature j was not selected then sel_mat[i, j] = 0.

    # n_sims: TODO: check if this argument is no longer needed?

    # calc_errors: logical; whether or not to calculate a confidence interval
    # for the stability metric.

    # conf: numeric; if calc_errors is TRUE, then a confidence interval will
    # be calculated at this confidence level. (Must be strictly between 0 and
    # 1.)

    # coarseness: TODO: check if this argument is no longer needed?

    # Returns:

    # If calc_errors is FALSE: returns the value
    # of the stability metric, which is at most 1. (Returns NA if metric
    # is unable to be calculated, for example if sel_mat contains all 0s.) If
    # calc_errors is TRUE, returns a length
    # 3 numeric vector containing the calculated stability metric, a lower
    # bound for the conf confidence interval, and an upper bound for the
    # confidence interval.

    if(is.na(n_sims)){
        n_sims <- nrow(sel_mat)/coarseness
    }

    # Check inputs
    checkNSBStabInputsNone(sel_mat, n_sims, calc_errors, conf, coarseness)
    p <- ncol(sel_mat)
    
    if(all(sel_mat == 0)){
        return(NA)
    }
    # Eliminate rows of all zeroes (no selections)
    all_zero_rows <- apply(sel_mat, 1, function(x){all(x == 0)})
    M <- sum(!all_zero_rows)
    if(M <= 1){
        return(NA)
    }

    if(M <= n_sims*coarseness/10){
        print("Warning here")
        warning("sum(!all_zero_rows) <= n_sims*coarseness/10")
    }
    sel_mat_nonzero <- sel_mat[!all_zero_rows, ]
    stopifnot(ncol(sel_mat_nonzero) == p)
    stopifnot(M == nrow(sel_mat_nonzero))

    # Calculate s_f_squared values
    p_hat <- colMeans(sel_mat_nonzero)
    s_f_squared <- apply(sel_mat_nonzero, 2,
        function(x){M/(M-1)*mean(x)*(1-mean(x))})
    stopifnot(all.equal(s_f_squared, M/(M-1)*p_hat*(1-p_hat)))
    k <- rowSums(sel_mat_nonzero)
    k_bar <- mean(k)

    stat <- 1 - mean(s_f_squared)/(k_bar/p*(1 - k_bar/p))

    # Check output
    stopifnot(length(stat) == 1)
    stopifnot(is.numeric(stat) | is.integer(stat))
    stopifnot(stat <= 1)

    if(calc_errors){
        alpha <- 1 - conf
        stopifnot(alpha >= 0 & alpha <= 1)

        # Calculate phi_hat_i's
        phi_hat_i <- numeric(M)
        stopifnot(ncol(sel_mat_nonzero) == length(p_hat))
        for(i in 1:M){
            phi_hat_i[i] <- 1/(k_bar/p*(1 - k_bar/p))*
                (mean(sel_mat_nonzero[i, ]*p_hat) - k[i]*k_bar/p^2 +
                stat/2*(2*k_bar*k[i]/p^2 - k[i]/p - k_bar/p + 1))
        }
        stopifnot(all(phi_hat_i != 0))
        v <- 4/M^2*sum((phi_hat_i - mean(phi_hat_i))^2)
        stopifnot(v >= 0)
        # Calculate errors for error bars
        margin <- qnorm(1 - alpha/2)*sqrt(v)

        # Check output

        stopifnot(length(margin) == 1)
        stopifnot(is.numeric(margin) | is.integer(margin))
        stopifnot(margin <= 1)
        stopifnot(margin >= 0)

        # Upper margin can't be greater than 1 since statistic can't be
        # greater than 1
        return(c(stat, stat - margin, min(stat + margin, 1)))
    }

    return(stat)

}

checkNSBStabInputs <- function(sel_mat, n_sims, calc_errors, conf,
    coarseness, cluster_count, n_clust_feats, C){
    stopifnot(is.matrix(sel_mat))
    stopifnot(all(!is.na(sel_mat)))
    stopifnot(all(sel_mat %in% c(0, 1)))

    stopifnot(length(n_sims) == 1)
    stopifnot(n_sims == round(n_sims))
    stopifnot(n_sims > 1)

    stopifnot(length(calc_errors) == 1)
    stopifnot(is.logical(calc_errors))

    stopifnot(length(conf) == 1)
    stopifnot(conf > 0)
    stopifnot(conf < 1)

    stopifnot(nrow(sel_mat) == n_sims*coarseness)

    stopifnot(length(cluster_count) == 1)
    stopifnot(is.character(cluster_count))
    if(cluster_count %in% c("single", "multi") & (n_clust_feats==0)){
        stop("Must provide n_clust_feats if cluster_count is single or multi.")
        stopifnot(length(n_clust_feats) == 1)
        stopifnot(n_clust_feats == round(n_clust_feats))
        stopifnot(n_clust_feats > 1)
    }

    stopifnot(cluster_count %in% c("none", "single", "multi", "corr"))

    p <- ncol(sel_mat)

    if(cluster_count == "corr"){
        if(any(is.na(C))){
            stop("Must provide C if cluster_count==corr")
        }
        stopifnot(ncol(C) == p)
        stopifnot(nrow(C) == p)
        stopifnot(all(C >= 0))
        stopifnot(all(C <= 1))
    }
}

checkNSBStabInputsNone <- function(sel_mat, n_sims, calc_errors, conf,
    coarseness){
    stopifnot(is.matrix(sel_mat))
    stopifnot(all(!is.na(sel_mat)))
    stopifnot(all(sel_mat %in% c(0, 1)))

    stopifnot(length(n_sims) == 1)
    stopifnot(n_sims == round(n_sims))
    stopifnot(n_sims > 1)

    stopifnot(length(calc_errors) == 1)
    stopifnot(is.logical(calc_errors))

    stopifnot(length(conf) == 1)
    stopifnot(conf > 0)
    stopifnot(conf < 1)

    stopifnot(nrow(sel_mat) == n_sims*coarseness)
}

createHamStabPlot <- function(stab_mets, n_model, p, k,
    beta_high, legend=TRUE){
    require(ggplot2)

    n_methods <- ncol(stab_mets)
    p_max <- nrow(stab_mets)
    j_choices <- 1:p_max
    names <- colnames(stab_mets)

    # Eliminate any rows with all NA entries
    all_na_rows <- apply(stab_mets, 1, function(x){all(is.na(x))})
    j_choices_to_eliminate <- which(all_na_rows)
    j_choices <- setdiff(j_choices, j_choices_to_eliminate)
    stab_mets <- stab_mets[!all_na_rows, ]

    df_gg <- data.frame(rep(j_choices, n_methods), as.vector(stab_mets),
        rep(names, each=length(j_choices)))

    colnames(df_gg) <- c("Model_Size", "HammingInstability", "Method")

    # Make labels display-friendly
    df_gg <- makeDfMethodsDisplay(df_gg)

    # Get max ranking value for integer labels on horizontal axis
    max_rank <- max(df_gg$Model_Size)
    if((max_rank %% 2) != 0){
        max_rank <- max_rank + 1
    }

    subtitle <- paste("n = ", n_model, ", p = ", p, ", k = ", k, ",
        beta_high = ", beta_high, sep="")

    plot <- ggplot(df_gg, aes(x=Model_Size, y=HammingInstability, color=Method,
        shape=Method)) + suppressWarnings(geom_point(size=2.5, alpha=1)) +
        xlab("Model Size") + ylab("Hamming Instability") +
        labs(subtitle=subtitle) +
        scale_x_continuous(breaks=seq(2, max_rank, by=2))

    if(!legend){
        plot <- plot + theme(legend.position="none")
    }

    return(plot)
}

getNSBDf <- function(stab_mets, plot_errors=FALSE, lowers=NA, uppers=NA){
    n_methods <- ncol(stab_mets)
    p_max <- nrow(stab_mets)
    j_choices <- 1:p_max
    names <- colnames(stab_mets)

    stopifnot(n_methods == length(names))

    # Eliminate any rows with all NA entries
    all_na_rows <- apply(stab_mets, 1, function(x){all(is.na(x))})
    j_choices_to_eliminate <- which(all_na_rows)
    j_choices <- setdiff(j_choices, j_choices_to_eliminate)
    stab_mets <- stab_mets[!all_na_rows, ]
    if(plot_errors){
        lowers <- lowers[!all_na_rows]
        uppers <- uppers[!all_na_rows]
    }

    stopifnot(length(as.vector(stab_mets)) == length(j_choices)*n_methods)

    df_gg <- data.frame(rep(j_choices, n_methods), as.vector(stab_mets),
        rep(names, each=length(j_choices)))

    colnames(df_gg) <- c("Model_Size", "NSBStability", "Method")

    # Make labels display-friendly
    df_gg <- makeDfMethodsDisplay(df_gg)

    if(plot_errors){
        # Add in error bounds
        lower_bounds <- as.vector(lowers)
        upper_bounds <- as.vector(uppers)
        stopifnot(length(lower_bounds) == nrow(df_gg) &
            length(upper_bounds) == nrow(df_gg))
        df_gg <- data.frame(df_gg, lower_bounds, upper_bounds)
        colnames(df_gg)[(ncol(df_gg) - 1):ncol(df_gg)] <- c("StabLower",
            "StabUpper")
    }

    return(df_gg)
}

createNSBStabPlot <- function(stab_mets, n_model=NA, p=NA, k=NA,
    beta_high=NA, legend=TRUE, plot_errors=FALSE, lowers=NA, uppers=NA,
    cluster_count="none", xaxis="No. Fitted Coefficients", subtitle=TRUE){
    require(ggplot2)

    stopifnot(cluster_count %in% c("none", "single", "multi", "corr"))

    if(subtitle & any(is.na(c(n_model, p, k, beta_high)))){
        stop("if subtitle is TRUE, must provide n_model, p, k, beta_high")
    }

    if(plot_errors){
        stopifnot(all.equal(dim(lowers), dim(uppers), dim(stab_mets)))
    }
    
    n_methods <- ncol(stab_mets)

    # Form data.frame
    df_gg <- getNSBDf(stab_mets, plot_errors, lowers, uppers)
    
    # Get max ranking value for integer labels on horizontal axis
    max_rank <- max(df_gg$Model_Size)
    if((max_rank %% 2) != 0){
        max_rank <- max_rank + 1
    }

    if(subtitle){
        subtitle_txt <- paste("n = ", n_model, ", p = ", p, ", k = ", k, ",
            beta_high = ", beta_high, sep="")
    }

    plot <- ggplot(df_gg, aes(x=Model_Size, y=NSBStability,
        color=Method, shape=Method)) + scale_shape_manual(values=1:n_methods)

    plot <- plot + suppressWarnings(geom_point(size=2.5, alpha=1)) +
        xlab(xaxis) + ylab("NSB Stability") +
        scale_x_continuous(breaks=seq(2, max_rank, by=2))

    if(subtitle){
        plot <- plot + labs(subtitle=subtitle_txt)
    }

    if(cluster_count=="multi"){
        plot <- plot + ylab("NSB Stability (cluster counting, multi)")
    } else if(cluster_count=="none"){
        plot <- plot + ylab("NSB Stability")
    } else if(cluster_count=="single"){
        plot <- plot + ylab("NSB Stability (cluster counting, single)")
    } else if(cluster_count=="corr"){
        plot <- plot + ylab("NSB Stability (2020 correlated metric)")
    }

    if(plot_errors){
        plot <- plot + geom_errorbar(aes(ymin = StabLower, ymax = StabUpper),
            width = 0.5)
    }

    if(!legend){
        plot <- plot + theme(legend.position="none")
    }

    return(plot)
}

createStabMSEPlot <- function(stab_mets, losses, names, j_choices_mat,
    n_model=NA, p=NA, k=NA, beta_high=NA, J=NA, legend=TRUE, plot_errors=FALSE,
    lowers=NA, uppers=NA, cluster_count="none", xaxis="No. Fitted Coefficients",
    cluster_count_mat=NA, subtitle=TRUE){

    require(ggplot2)
    require(dplyr)

    stopifnot(cluster_count %in% c("none", "single", "multi", "corr"))

    if(subtitle & any(is.na(c(n_model, p, k, beta_high)))){
        stop("if subtitle is TRUE, must provide n_model, p, k, beta_high")
    }

    n_methods <- ncol(stab_mets)

    stopifnot(n_methods == length(names))

    n_sims <- ncol(losses)

    # If J is not provided, figure out on own
    stopifnot(length(J) <= 1)
    if(is.na(J)){
        J <- nrow(losses)/n_methods
    }
    stopifnot(J == round(J))
    stopifnot(!(any(losses[!is.na(losses)] < 0)))
    stopifnot(!(n_sims != ncol(j_choices_mat) | nrow(losses) != nrow(j_choices_mat)))

    if(plot_errors){
        stopifnot(all.equal(dim(lowers), dim(uppers), dim(stab_mets)))
    }

    if(xaxis=="No. Selected Features"){
        stopifnot(!(n_sims != ncol(cluster_count_mat) | nrow(losses) != nrow(cluster_count_mat)))
        stopifnot(all(cluster_count_mat[j_choices_mat] != 0))
    }

    if(xaxis=="No. Selected Features" & any(is.na(cluster_count_mat))){
        stop("If xaxis=No. Selected Features, must provide cluster_count_mat")
    }

    # Form data.frame
    if(cluster_count == "none"){
        df_gg_stab <- getNSBDf(stab_mets, plot_errors, lowers, uppers)
    } else{
        df_gg_stab <- getNSBDf(stab_mets, plot_errors=FALSE, lowers=NA,
            uppers=NA)
    }
    
    df_gg_mse <- getErrorDf(losses, names, j_choices_mat, J, xaxis,
        cluster_count_mat, plot_errors=plot_errors)

    df_gg <- left_join(df_gg_stab, df_gg_mse, by=c("Model_Size", "Method"))

    if(subtitle){
        subtitle_txt <- paste("n = ", n_model, ", p = ", p, ", k = ", k, ",
            beta_high = ", beta_high, sep="")
    }

    plot <- ggplot(df_gg, aes(x=NSBStability, y=MSE,
        color=Method, shape=Method)) + scale_shape_manual(values=1:n_methods) +
        suppressWarnings(geom_point(size=2.5, alpha=1)) 
        # + 
        # geom_hline(yintercept=true_model_loss
        #     # , color="red"
        #     , linetype = "dashed"
        #     )

    if(subtitle){
        plot <- plot + labs(subtitle=subtitle_txt)
    }

    if(cluster_count=="multi"){
        plot <- plot + xlab("NSB Stability (cluster counting, multi)")
    } else if(cluster_count=="none"){
        plot <- plot + xlab("NSB Stability")
    } else if(cluster_count=="single"){
        plot <- plot + xlab("NSB Stability (cluster counting, single)")
    } else if(cluster_count=="corr"){
        plot <- plot + xlab("NSB Stability (2020 correlated metric)")
    }

    if(plot_errors){
        plot <- plot + geom_errorbar(aes(ymin = MSELower, ymax = MSEUpper)
            , width = 0.05
            )
        if(cluster_count == "none"){
            plot <- plot + geom_errorbar(aes(xmin = StabLower, xmax = StabUpper)
            , width = 0.25
            )
        }
    }

    if(!legend){
        plot <- plot + theme(legend.position="none")
    }

    return(plot)
}


saveFigure <- function(subdir, plot, size="large", filename=NA){

    stopifnot(length(filename) == 1)

    stopifnot(size %in% c("small", "medium", "mlarge", "large", "xlarge",
        "slide"))

    dir_fig <- file.path(folder_dir, subdir)
    dir.create(dir_fig, showWarnings = FALSE, recursive = TRUE)
    setwd(dir_fig)

    if(is.na(filename)){
        filename <- paste("p=", p_p_prime, ",beta=", round(beta_high_p_prime, 2),
            sep="")
        filename <- gsub(pattern="\\.", replacement="_", x=filename)
        filename <- paste(filename, ".pdf", sep="")
    }

    if(size=="large"){
        pdf(file=filename, title=filename, paper="special", width=10.5,
            height=5.5, bg="white")
    } else if(size=="mlarge"){
        pdf(file=filename, title=filename, paper="special", width=10.5,
            height=4, bg="white")
    } else if(size=="xlarge"){
        pdf(file=filename, title=filename, paper="special", width=10.5,
            height=6.5, bg="white")
    } else if(size=="small"){
        pdf(file=filename, title=filename, paper="special", width=5, height=5,
            bg="white")
    } else if(size=="medium"){
        pdf(file=filename, title=filename, paper="special", width=6.4,
            height=5.5, bg="white")
    } else if(size=="slide"){
        pdf(file=filename, title=filename, paper="special", width=540/72,
            height=360/72, bg="white")
    }
    
    plot(plot)
    dev.off()

    setwd(sim_dir)
}

get_legend<-function(myggplot){
    # From http://www.sthda.com/english/wiki/wiki.php?id_contents=7930#add-a-common-legend-for-multiple-ggplot2-graphs
    require(gridExtra)
    tmp <- ggplot_gtable(ggplot_build(myggplot))
    leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
    legend <- tmp$grobs[[leg]]
    return(legend)
}

partial_out <- function(X, y){
    # Regresses y against the features in X and returns the residual.

    # Check inputs
    stopifnot(is.matrix(X) | is.numeric(X))

    if(is.matrix(X)){
        n <- nrow(X)
        p <- ncol(X)
    } else{
        n <- length(X)
        p <- 1
    }

    stopifnot(n == length(y))
    stopifnot(is.numeric(y) | is.integer(y))

    if(is.matrix(X)){
        y_reg_df <- as.data.frame(cbind(y, X))
        colnames(y_reg_df)[1] <- "y"
    } else{
        y_reg_df <- data.frame("y" = y, X)
    }
    
    y_reg <- lm(y ~ ., y_reg_df)
    y_hat <- predict(y_reg)
    stopifnot(is.numeric(y_hat))
    y_residual <- y - y_hat

    # Check output

    stopifnot(length(y_residual) == n)
    stopifnot(is.numeric(y_residual))

    return(y_residual)
}





# TODO: After simplifying discardClusterMembers, change this function to that
# only takes in colMeans(res), which is all that is needed in
# discardClusterMembers
getSelected <- function(res, clus_sel_props_p, clusters, cutoff=0){
    # Returns p-vector containing indices of features of x sorted in 
    # decreasing order of selection proportions, based on matrices of selection
    # indcators. Also returns clusters corresponding to features and whether or
    # not a given feature is tied in selection proportion with the next feature.

    # Inputs:

    # res
    # a B x p (or 2*B x p for "SS") numeric matrix containing indicator
    # variables; res[i, j] = 1 if feature j was selected on subsample i, and 0
    # otherwise.

    # clus_sel_props_p
    # a p-vector; entry j is the selection proportion for the cluster
    # corresponding to feature j.

    # clusters
    # A list of vectors of integers representing clusters of similar features
    # (derived by function formatClusters from matrix R that is input to
    # css).

    # cutoff
    # Numeric; optional argument. If provided, css will return only those
    # clusters with selection proportions equal to at least cutoff. Must be
    # between 0 and 1.

    # Outputs

    # selected
    # an integer vector of length at most p containing indices of features of x 
    # sorted in decreasing order of selection proportions. (Length may be less
    # than p because of clusters; only one member of each cluster will appear
    # in selected.)

    # Check inputs

    stopifnot(is.matrix(res))
    stopifnot(all(!is.na(res)))
    stopifnot(all(res %in% c(0, 1)))

    p <- ncol(res)

    stopifnot(is.numeric(clus_sel_props_p))
    stopifnot(length(clus_sel_props_p) == p)
    stopifnot(all(clus_sel_props_p >= 0))
    stopifnot(all(clus_sel_props_p <= 1))
    stopifnot(p == ncol(res))

    stopifnot(all(lengths(clusters) >= 1))
    stopifnot(is.list(clusters))
    stopifnot(all(!is.na(clusters)))
    stopifnot(length(clusters) == length(unique(clusters)))

    if(length(clusters) >= 1){
        for(i in 1:length(clusters)){
            stopifnot(length(clusters[[i]]) >= 1)
            stopifnot(length(clusters[[i]]) == length(unique(clusters[[i]])))
            stopifnot(all(!is.na(clusters[[i]])))
            stopifnot(is.integer(clusters[[i]]))
        }

        if(length(clusters) >= 2){
            # Check that clusters are non-overlapping
            for(i in 1:(length(clusters) - 1)){
                for(j in (i+1):length(clusters)){
                    stopifnot(length(intersect(clusters[[i]],
                        clusters[[j]])) == 0)
                }
            }
        }
    }

    stopifnot(is.numeric(cutoff) | is.integer(cutoff))
    stopifnot(length(cutoff) == 1)
    stopifnot(cutoff >= 0)
    stopifnot(cutoff <= 1)
    
    # Sort selected in descending order of selection proportions
    selected <- (1:p)[order(clus_sel_props_p, decreasing=TRUE)]

    # If cutoff was specified, remove features whose cluster selection
    # proportions are below the cutoff
    if(cutoff > 0){
        clus_sel_props_p_sorted <- sort(clus_sel_props_p, decreasing=TRUE)
        selected <- selected[clus_sel_props_p_sorted >= cutoff]
        # Possible that selected is now empty; if so, return an empty set.
        if(length(selected) == 0){
            return(list(selected=selected, clusters=clusters))
        }
    }

    # selected now contains all p features (or all features corresponding to
    # clusters with selection proportions greater than cutoff, if cutoff was
    # specified), but we only want one entry 
    # corresponding to each cluster. In order to do that, we will remove all
    # but one of the features from each cluster from selected. (We will choose
    # the cluster member that was selected
    # the most frequently individually, which we can discern from res.)

    if(length(clusters) > 0){
        stopifnot(all(lengths(clusters) >= 1))
        for(i in 1:length(clusters)){
            if(length(clusters[[i]]) > 1){
                selected <- discardClusterMembers(clusters[[i]], selected,
                    res)
            }  
        }
    }

    # Check outputs

    final_sel_props <- clus_sel_props_p[selected]

    p_selected <- length(selected)

    stopifnot(p_selected != 0)
    stopifnot(length(final_sel_props) == p_selected)

    for(j in 1:(p_selected - 1)){
        stopifnot(final_sel_props[j] >= final_sel_props[j + 1])
    }

    stopifnot(is.integer(selected))
    stopifnot(length(unique(selected)) == p_selected)
    stopifnot(p_selected <= p)

    ret <- list(selected=selected)

    return(ret)
}

# TODO: change this function to that only takes in colMeans(res)
discardClusterMembers <- function(cluster_i, selected, res){
    # Given a set of selected features selected and a particular cluster
    # cluster_i, removes all but the most frequently selected cluster member
    # in cluster_i from selected.

    stopifnot(is.integer(cluster_i))
    stopifnot(length(cluster_i) > 1)
    stopifnot(!is.list(cluster_i))

    most_sel <- integer()

    if(any(cluster_i %in% selected)){
        # Identify the most frequently selected cluster member
        sel_props <- colMeans(res)[cluster_i]
        most_sel <- cluster_i[which.max(sel_props)]
        stopifnot(length(most_sel) == 1)
        # Remove all but the most frequently selected clustered feature 
        # from selected
        selected <- setdiff(selected, setdiff(cluster_i,
            most_sel))
        stopifnot(most_sel %in% selected)
    }

    # Check output
    stopifnot(!(setdiff(cluster_i, most_sel) %in% selected))

    return(as.integer(selected))
}

selAveraging <- function(selected, clusters, feat_sel_props, weighting){
    # Calculates weights for each cluster member. In particular, for each
    # selected feature, this function 
    # provides (1) a logical as to whether or not this feature corresponds to a 
    # cluster and requires averaging, (2) a vector of indices of the features in 
    # the cluster containing that feature, and (3) a vector of weights to apply 
    # to the features.
    # Inputs

    # selected
    # an integer vector of length at most p containing indices of features of x 
    # sorted in decreasing order of selection proportions. (Length may be less
    # than p because of clusters; only one member of each cluster will appear
    # in selected.)

    # clusters
    # A list of vectors of integers representing clusters of similar 
    # features (derived from R). Only contains clusters of size 
    # greater than 1.

    # res
    # a B x p (or 2*B x p for "SS") integer matrix. Each row of res is a 
    # selected set (res_ij is 1 if feature j was selected on resample i 
    # and 0 if not).

    # weighting
    # character specifying weighting scheme



    # Outputs:

    # to_avg: A logical vector to_avg of the same length as 
    # selected. If feature j is in a cluster, the jth entry of to_avg will be
    # TRUE; otherwise will be FALSE. 
    
    # avg_feats: A list of the same length as selected of integer vectors. 
    # avg_feats[[j]] is a vector of indices of the features from the
    # cluster of which feature j is a member.
    
    # weights: A list of the same length as selected of numeric vectors.
    # weights[j] is the weights to use (in the same order as the jth entry of
    # avg_feats).

    # Check inputs

    stopifnot(is.integer(selected))

    p_ret <- length(selected)

    stopifnot(is.list(clusters))
    stopifnot(all(lengths(clusters) >= 1))
    stopifnot(all(!is.na(clusters)))
    stopifnot(length(clusters) == length(unique(clusters)))
    stopifnot(length(clusters) == p_ret)

    if(length(clusters) >= 1){
        for(i in 1:length(clusters)){
            stopifnot(length(clusters[[i]]) == length(unique(clusters[[i]])))
            stopifnot(all(!is.na(clusters[[i]])))
            stopifnot(is.integer(clusters[[i]]))
        }
        if(length(clusters) >= 2){
            # Check that clusters are non-overlapping
            for(i in 1:(length(clusters) - 1)){
                for(j in (i+1):length(clusters)){
                    stopifnot(length(intersect(clusters[[i]],
                        clusters[[j]])) == 0)
                }
            }
        }
    }

    stopifnot(length(weighting)==1)
    stopifnot(is.character(weighting))
    stopifnot(weighting %in% c("sparse", "simple_avg", "weighted_avg"))


    # Identify weights
    to_avg <- logical(p_ret)
    avg_feats <- list()
    weights <- list()
    if(p_ret > 0){
        for(j in 1:p_ret){
            # Find the cluster feature j is a member of
            j_clust_found <- FALSE
            for(i in 1:length(clusters)){
                cluster_i <- clusters[[i]]
                name_i <- names(clusters)[i]
                stopifnot(is.integer(cluster_i))

                # Find the cluster this feature is in (if it's not cluster_i,
                # move on to the next cluster)
                if(selected[j] %in% cluster_i){
                    results <- getWeights(cluster_i, name_i, j, to_avg,
                        avg_feats, weights, weighting, feat_sel_props)

                    to_avg <- results$to_avg
                    avg_feats <- results$avg_feats
                    weights <- results$weights

                    rm(results)

                    # End this for loop (don't need to check remaining
                    # clusters, since we already found the cluster this
                    # feature is in)
                    j_clust_found <- TRUE
                    break
                }
            }
            if(!j_clust_found){
                stop(paste("No cluster found for feature", selected[j]))
            }
            # # This code is only reached if feature j is not in a cluster, in
            # # which case it is assumed to be in a "cluster" containing only
            # # itself.
            # results <- getWeights(selected[j],
            #     paste("c", selected[j], sep=""), j, to_avg, avg_feats,
            #     weights, weighted, res)

            # to_avg <- results$to_avg
            # avg_feats <- results$avg_feats
            # weights <- results$weights

            # rm(results)
        }
    }

    # Check output

    stopifnot(length(to_avg) == p_ret)
    stopifnot(is.logical(to_avg))
    stopifnot(all(!is.na(to_avg)))
    stopifnot(all(to_avg))

    stopifnot(is.list(avg_feats))
    stopifnot(length(avg_feats) <= p_ret)

    stopifnot(length(weights) == p_ret)
    stopifnot(p_ret == length(avg_feats))

    for(i in 1:p_ret){
        if(to_avg[i]){
            stopifnot(is.integer(avg_feats[[i]]))
            stopifnot(length(avg_feats[[i]]) >= 1)
            stopifnot(all(!is.na(avg_feats[[i]])))
            stopifnot(length(unique(avg_feats[[i]])) ==
                    length(avg_feats[[i]]))
            stopifnot(selected[i] %in% avg_feats[[i]])

            stopifnot(length(avg_feats[[i]]) == length(weights[[i]]))
            stopifnot(all(weights[[i]] >= 0))
            stopifnot(all(weights[[i]] <= 1))
            stopifnot(abs(sum(weights[[i]]) - 1) < 10^(-6))
        } else{
            # If there are any TRUE values later in to_avg, then the ith
            # entries of weights and avg_feats should be empty. If not,
            # then weights and avg_feats should have lengths shorter than i.
            if(any(to_avg[(i + 1):p_ret])){
                stopifnot(length(avg_feats[[i]]) == 0)
                stopifnot(length(weights[[i]]) == 0)
            } else{
                stopifnot(i > p_ret)
            }
        }
    }

    stopifnot(is.list(avg_feats))
    stopifnot(length(avg_feats) == length(selected))
    stopifnot(!is.null(names(avg_feats)))
    stopifnot(all(!is.na(avg_feats) & (avg_feats != "")))

    stopifnot(is.list(weights))
    stopifnot(length(avg_feats) == length(weights))

    if(length(avg_feats) > 0){
        for(i in 1:length(avg_feats)){
            stopifnot(is.integer(avg_feats[[i]]))
            stopifnot(length(avg_feats[[i]]) >= 1)
            stopifnot(all(!is.na(avg_feats[[i]])))
            stopifnot(length(unique(avg_feats[[i]])) ==
                length(avg_feats[[i]]))
            stopifnot(selected[i] %in% avg_feats[[i]])

            stopifnot(length(avg_feats[[i]]) == length(weights[[i]]))
            stopifnot(is.numeric(weights[[i]]))
            stopifnot(all(weights[[i]] >= 0))
            stopifnot(all(weights[[i]] <= 1))
            stopifnot(abs(sum(weights[[i]]) - 1) < 10^(-6))
        }
    }

    # Add names to weights
    names(weights) <- names(avg_feats)

    return(list(to_avg=to_avg, avg_feats=avg_feats, weights=weights))
}

getWeights <- function(cluster_i, name_i, j, to_avg, avg_feats, weights,
    weighting, feat_sel_props){

    to_avg[j] <- TRUE
    avg_feats[[j]] <- cluster_i
    names(avg_feats)[j] <- name_i

    # Get the selection proportions of each cluster member
    sel_props <- feat_sel_props[cluster_i]
    stopifnot(all(sel_props >= 0))
    stopifnot(all(sel_props <= 1))

    n_weights <- length(cluster_i)
    
    # Weighted or simple average?
    if(weighting == "sparse"){
        # Sparse cluster stability selection: All features in cluster with
        # selection proportion equal to the max
        # for the cluster get equal weight; rest of cluster gets 0 weight
        if(sum(sel_props) == 0){
            weights_i <- rep(1/n_weights, n_weights)
        } else{
            maxes <- sel_props==max(sel_props)
            stopifnot(sum(maxes) > 0)
            stopifnot(sum(maxes) <= n_weights)
            weights_i <- rep(0, n_weights)
            weights_i[maxes] <- 1/sum(maxes)
        }
    } else if(weighting == "weighted_avg"){
        # Get weights for weighted average
        if(sum(sel_props) == 0){
            weights_i <- rep(1/n_weights, n_weights)
        } else{
            weights_i <- sel_props/sum(sel_props)
        }
    } else if(weighting == "simple_avg"){
        weights_i <- rep(1/n_weights, n_weights)
    } else{
        stop("weighting must be one of sparse, simple_avg, or weighted_avg")
    }

    stopifnot(abs(sum(weights_i) - 1) < 10^(-6))
    stopifnot(length(weights_i) == length(cluster_i))
    stopifnot(length(weights_i) >= 1)

    weights[[j]] <- weights_i

    return(list(to_avg=to_avg, avg_feats=avg_feats, weights=weights))
}





cssWrapper <- function(X, y, lambda
    , clusters = list()
    , cutoff = 0
    , fitfun = cssLasso
    , sampling_type = "SS"
    , B = ifelse(sampling_type == "MB", 100L, 50L)
    , weighting="sparse"
    # , average = FALSE
    # , weighted = TRUE
    , prop_feats_remove = 0
    # , ...
    ){

    # Check inputs

    feat_names <- as.character(NA)
    if(!is.null(colnames(X))){
        feat_names <- colnames(X)
    }
    clust_names <- as.character(NA)
    if(!is.null(names(clusters))){
        clust_names <- names(clusters)
    }

    # Check if x is a matrix; if it's a data.frame, convert to matrix.
    if(is.data.frame(X)){
        X <- model.matrix(~ ., X)
    }

    stopifnot(is.matrix(X))
    stopifnot(all(!is.na(X)))

    n <- nrow(X)
    p <- ncol(X)
    stopifnot(p >= 2)
    if(length(feat_names) > 1){
        stopifnot(length(feat_names) == p)
    } else{
        stopifnot(is.na(feat_names))
    }
    
    colnames(X) <- character()

    stopifnot(length(y) == n)
    # Intentionally don't check y or lambda further to allow for flexbility--these
    # inputs should be checked within fitfun.

    stopifnot(!is.na(clusters))
    if(is.list(clusters)){
        stopifnot(all(!is.na(clusters)))
        stopifnot(length(clusters) == length(unique(clusters)))

        if(is.list(clusters) & length(clusters) > 0){
            for(i in 1:length(clusters)){
                stopifnot(length(clusters[[i]]) == length(unique(clusters[[i]])))
                stopifnot(all(!is.na(clusters[[i]])))
                stopifnot(is.integer(clusters[[i]]) | is.numeric(clusters[[i]]))
                stopifnot(all(clusters[[i]] == round(clusters[[i]])))
                clusters[[i]] <- as.integer(clusters[[i]])
            }

            if(length(clusters) >= 2){
                # Check that clusters are non-overlapping
                for(i in 1:(length(clusters) - 1)){
                    for(j in (i+1):length(clusters)){
                        if(length(intersect(clusters[[i]], clusters[[j]])) != 0){
                            error_mes <- paste("Overlapping clusters detected; clusters must be non-overlapping. Overlapping clusters: ",
                                i, ", ", j, ".", sep="")
                            stop(error_mes)
                        }
                    }
                }
            }
        } 
    } else{
        stopifnot(is.numeric(clusters) | is.integer(clusters))
        stopifnot(length(clusters) == length(unique(clusters)))
        stopifnot(all(!is.na(clusters)))
        stopifnot(is.integer(clusters) | is.numeric(clusters))
        stopifnot(all(clusters == round(clusters)))
        clusters <- as.integer(clusters)

    }

    stopifnot(is.numeric(cutoff) | is.integer(cutoff))
    stopifnot(length(cutoff) == 1)
    stopifnot(cutoff >= 0)
    stopifnot(cutoff <= 1)

    stopifnot(class(fitfun) == "function")
    stopifnot(length(fitfun) == 1)
    if(!identical(formals(fitfun), formals(cssLasso))){
        err_mess <- paste("fitfun must accept arguments named X, y, and lambda. Detected arguments to fitfun:",
            paste(names(formals(fitfun)), collapse=", "))
        stop(err_mess)
    }
    
    stopifnot(is.character(sampling_type))
    stopifnot(length(sampling_type) == 1)
    stopifnot(sampling_type %in% c("SS", "MB"))

    stopifnot(length(B) == 1)
    stopifnot(is.numeric(B) | is.integer(B))
    stopifnot(B == round(B))
    stopifnot(B > 0)
    if(B < 10){
        warning("Small values of B may lead to poor results.")
    } else if (B > 2000){
        warning("Large values of B may require long computation times.")
    }

    stopifnot(length(weighting)==1)
    if(!is.character(weighting)){
        stop("Weighting must be a character")
    }
    if(!(weighting %in% c("sparse", "simple_avg", "weighted_avg"))){
        stop("Weighting must be a character and one of sparse, simple_avg, or weighted_avg")
    }

    stopifnot(length(prop_feats_remove) == 1)
    stopifnot(is.numeric(prop_feats_remove) | is.integer(prop_feats_remove))
    stopifnot(prop_feats_remove >= 0 & prop_feats_remove < 1)
    if(prop_feats_remove > 0){
        # Make sure p is at least 2 or else this doesn't make sense
        stopifnot(p >= 2)
    }

    ### Create subsamples

    subsamps_object <- createSubsamples(n, p, B, sampling_type,
        prop_feats_remove)

    ### Get matrix of selected feature sets from subsamples

    res <- getSelMatrix(X, y, lambda, B, sampling_type, subsamps_object, fitfun)

    stopifnot(is.matrix(res))
    if(sampling_type=="SS"){
        stopifnot(nrow(res) == 2*B)
    } else{
        stopifnot(nrow(res) == B)
    }
    stopifnot(ncol(res) == p)
    stopifnot(all(res %in% c(0, 1)))

    #############################################

    # Format clusters into a list with no clusters of size 1
    cluster_res <- formatClusters(clusters, p=p, clust_names=clust_names)

    clusters <- cluster_res$clusters
    # clust_to_name_dict <- cluster_res$clust_to_name_dict

    rm(cluster_res)

    ### Get selection proportions for clusters corresponding to each feature

    clus_prop_results <- getClusterProps(clusters, res, sampling_type
        # , clust_to_name_dict
        )

    feat_sel_props <- clus_prop_results$feat_sel_props
    clus_sel_props_p <- clus_prop_results$clus_sel_props_p
    res_clus_p <- clus_prop_results$res_clus_p
    res_n_clusters <- clus_prop_results$res_n_clusters

    rm(clus_prop_results)

    # Get selected set (sort features in decreasing order of weighted selection
    # proportions, including only one feature per cluster)

    selected_results <- getSelected(res, clus_sel_props_p, clusters, cutoff)

    selected <- selected_results$selected
    # multiple <- selected_results$multiple
    # clusters <- selected_results$clusters

    rm(selected_results)

    # If we are using averaging, in addition to the ordered selected set,
    # will also return lists avg_feats and weights of the same length.
    # If feature j is in a cluster, the jth entry of avg_feats will contain
    # the features from the cluster
    # of which feature j is a member, and the jth entry of weights will
    # be the weights to use (in the same order as the jth entry of
    # avg_feats).
    sel_avg_res <- selAveraging(selected, clusters, feat_sel_props, weighting)

    selected_clusts <- sel_avg_res$avg_feats
    weights <- sel_avg_res$weights

    rm(sel_avg_res)


    # Check outputs

    stopifnot(is.integer(selected))
    stopifnot(length(selected) <= p)
    stopifnot(length(selected) == length(unique(selected)))
    if(length(selected) > 0){
        stopifnot(max(selected) <= p)
        stopifnot(min(selected) >= 1)
    }

    stopifnot(is.numeric(feat_sel_props))
    stopifnot(length(feat_sel_props) == p)
    stopifnot(all(feat_sel_props >= 0))
    stopifnot(all(feat_sel_props <= 1))
    stopifnot(identical(feat_sel_props, colMeans(res)))

    stopifnot(length(B) == 1)
    stopifnot(is.numeric(B) | is.integer(B))
    stopifnot(B == round(B))
    stopifnot(B > 0)

    if(any(!is.na(feat_names))){
        names(selected) <- feat_names[selected]
        names(feat_sel_props) <- feat_names
        colnames(res) <- feat_names
        colnames(X) <- feat_names
    }

    stopifnot(ncol(res_n_clusters) == length(selected))


    stopifnot(is.list(selected_clusts))
    stopifnot(length(selected_clusts) == length(selected))
    stopifnot(!is.null(names(selected_clusts)))
    stopifnot(all(!is.na(selected_clusts) & (selected_clusts != "")))

    stopifnot(is.list(weights))
    stopifnot(length(selected_clusts) == length(weights))

    if(length(selected_clusts) > 0){
        for(i in 1:length(selected_clusts)){
            stopifnot(is.integer(selected_clusts[[i]]))
            stopifnot(length(selected_clusts[[i]]) >= 1)
            stopifnot(all(!is.na(selected_clusts[[i]])))
            stopifnot(length(unique(selected_clusts[[i]])) ==
                length(selected_clusts[[i]]))
            stopifnot(selected[i] %in% selected_clusts[[i]])

            stopifnot(length(selected_clusts[[i]]) == length(weights[[i]]))
            stopifnot(is.numeric(weights[[i]]))
            stopifnot(all(weights[[i]] >= 0))
            stopifnot(all(weights[[i]] <= 1))
            stopifnot(abs(sum(weights[[i]]) - 1) < 10^(-6))
        }
    }

    # Add names to weights
    names(weights) <- names(selected_clusts)

    # if(!is.na(clust_names)){
    #     # Add cluster names to selected_clusts, and weights
    #     names(selected_clusts) <- getClustNames(selected_clusts, clusters,)
    # }

    ret <- list(selected = selected,
        selected_clusts = selected_clusts,
        weights = weights,
        feat_sel_props = feat_sel_props,
        feat_sel_mat = res,
        clus_sel_props = colMeans(res_n_clusters),
        clus_sel_mat = res_n_clusters
        )


    class(ret) <- "cssr"
    
    return(ret)
}















    # TODO: figure out how to eliminate need for prototypes argument in
    # getClusterSelsFromGlmnet under clusterRepLasso (shouldn't be necessary for
    # anything)
getClusterSelsFromGlmnet <- function(lasso_sets, clusters, prototypes,
    n_cluster_members, non_cluster_feats, var_names_provided=FALSE,
    var_names=NA, averaging=FALSE){

    # Check inputs

    stopifnot(!is.list(clusters) | all(lengths(clusters) >= 1))
    stopifnot(is.list(clusters) | length(clusters) >= 1)

    stopifnot(is.integer(prototypes))
    if(is.list(clusters)){
        stopifnot(length(prototypes) == length(clusters))
    } else{
        stopifnot(length(prototypes) == 1)
    }
    stopifnot(all(!is.na(prototypes)))
    stopifnot(length(prototypes) == length(unique(prototypes)))

    stopifnot(is.numeric(n_cluster_members) | is.integer(n_cluster_members))
    stopifnot(n_cluster_members == round(n_cluster_members))
    stopifnot(n_cluster_members >= 0)
    # stopifnot(n_cluster_members <= p)

    stopifnot(is.numeric(non_cluster_feats) | is.integer(non_cluster_feats))
    # stopifnot(length(non_cluster_feats) <= p)
    stopifnot(length(non_cluster_feats) >= 0)
    if(length(non_cluster_feats) >= 1){
        stopifnot(length(non_cluster_feats) == length(unique(non_cluster_feats)))
        stopifnot(all(round(non_cluster_feats) == non_cluster_feats))
        # stopifnot(all(non_cluster_feats) %in% 1:p)
    }

    stopifnot(length(var_names_provided) == 1)
    stopifnot(is.logical(var_names_provided))
    if(var_names_provided){
        stopifnot(is.character(var_names))
        if(any(is.na(var_names))){
            stop("must provide var_names (with no NAs) if var_names_provided=TRUE")
        }
    }

    stopifnot(length(averaging) == 1)
    stopifnot(is.logical(averaging))

    if(is.list(clusters)){
        n_clusters <- length(clusters)
    } else{
        n_clusters <- 1
    }
    
    stopifnot(length(prototypes) == n_clusters)




    max_length <- max(vapply(lasso_sets, length, integer(1)))

    selected_sets <- list()
    if(averaging){
        to_avg_list <- list()
        selected_clusts_list <- list()
        weights_list <- list()
    }
    

    for(j in 1:max_length){
        # Lasso selected set of size j
        lasso_sets_j <- lasso_sets[lapply(lasso_sets, length) == j]
        if(length(lasso_sets_j) > 0){
            lasso_set_j <- lasso_sets_j[[1]]
            stopifnot(length(unique(lasso_set_j)) == length(lasso_set_j))

            # Recover features from original feature space: deal with
            # prototype and non-prototype features separately
            proto_inds <- lasso_set_j %in% 1:n_clusters
            lasso_set_j_protos <- lasso_set_j[proto_inds]
            lasso_set_j_non_protos <- lasso_set_j[!proto_inds]

            # Recover indices of non-prototype features

            if(length(lasso_set_j_non_protos) > 0){
                stopifnot(length(lasso_set_j_non_protos) <=
                    length(non_cluster_feats))
                stopifnot(all(lasso_set_j_non_protos - n_clusters %in%
                    1:length(non_cluster_feats)))

                recovered_feats <- non_cluster_feats[lasso_set_j_non_protos -
                    n_clusters]

                stopifnot(length(lasso_set_j_non_protos) ==
                    length(recovered_feats))
                stopifnot(all.equal(!proto_inds, lasso_set_j %in%
                    lasso_set_j_non_protos))
                stopifnot(length(recovered_feats) ==
                    length(unique(recovered_feats)))

                lasso_set_j_non_protos <- recovered_feats

                stopifnot(length(lasso_set_j_non_protos) == sum(!proto_inds))
            }

            stopifnot(length(unique(lasso_set_j_protos)) == length(lasso_set_j_protos))
            stopifnot(length(lasso_set_j_protos) == sum(proto_inds))

            proto_selected_j <- length(lasso_set_j_protos) > 0

            stopifnot(all(lasso_set_j_protos %in% 1:n_clusters))

            if(proto_selected_j){
                lasso_set_j_protos <- prototypes[lasso_set_j_protos]
                stopifnot(length(lasso_set_j_protos) == sum(proto_inds))
            }

            stopifnot(length(lasso_set_j) == length(lasso_set_j_protos) +
                length(lasso_set_j_non_protos))

            lasso_set_j[proto_inds] <- lasso_set_j_protos
            lasso_set_j[!proto_inds] <- lasso_set_j_non_protos

            stopifnot(length(unique(lasso_set_j)) == length(lasso_set_j))
            
            if(var_names_provided){
                stopifnot(max(lasso_set_j) <= length(var_names))
                selected_sets[[j]] <- var_names[lasso_set_j]
            } else{
                selected_sets[[j]] <- lasso_set_j
            }

            if(averaging){
                to_avg_list[[j]] <- logical(j)
                selected_clusts_list[[j]] <- list()
                weights_list[[j]] <- list()
            
                if(proto_selected_j){
                    ind_j <- matrix(FALSE, nrow=length(lasso_set_j),
                        ncol=n_clusters)
                    if(is.list(clusters)){
                        for(k in 1:n_clusters){
                            n_cluster_members_k <- length(clusters[[k]])

                            stopifnot(n_cluster_members_k >= 1)

                            ind_j[, k] <- lasso_set_j == prototypes[k]

                            stopifnot(sum(ind_j[, k]) %in% c(0, 1))

                            if(sum(ind_j[, k]) == 1){
                                to_avg_list[[j]][which(ind_j[, k])] <- TRUE

                                # selected_sets[[j]][lasso_set_j == 1] <- cluster_prototype
                                if(var_names_provided){
                                    selected_clusts_list[[j]][[which(ind_j[, k])]] <-
                                        var_names[clusters[[k]]]
                                } else{
                                    selected_clusts_list[[j]][[which(ind_j[, k])]] <-
                                        clusters[[k]]
                                }

                                weights_list[[j]][[which(ind_j[, k])]] <-
                                    rep(1/n_cluster_members_k, n_cluster_members_k)
                            }
                        }
                    } else{
                        stopifnot(n_clusters <= 1)
                        n_cluster_members_k <- length(clusters)
                        stopifnot(n_cluster_members_k != 1)

                        ind_j[, 1] <- lasso_set_j == prototypes

                        stopifnot(sum(ind_j[, 1]) %in% c(0, 1))

                        if(sum(ind_j[, 1]) == 1){
                            to_avg_list[[j]][which(ind_j[, 1])] <- TRUE

                            # selected_sets[[j]][lasso_set_j == 1] <- cluster_prototype
                            if(var_names_provided){
                                selected_clusts_list[[j]][[which(ind_j[, 1])]] <-
                                    var_names[clusters]
                            } else{
                                selected_clusts_list[[j]][[which(ind_j[, 1])]] <-
                                    clusters
                            }


                            
                            weights_list[[j]][[which(ind_j[, 1])]] <-
                                rep(1/n_cluster_members_k, n_cluster_members_k)

                            stopifnot(length(weights_list[[j]][[which(ind_j[, 1])]]) ==
                                length(selected_clusts_list[[j]][[which(ind_j[, 1])]]))
                        }
                    }
                    
                    stopifnot(sum(ind_j) != 0)
                    stopifnot(all(rowSums(ind_j) <= 1))
                } 
            }
        }
    }

    if(averaging){
        return(list(selected_sets=selected_sets,
            to_avg_list=to_avg_list, selected_clusts_list=selected_clusts_list,
            weights_list=weights_list))
    }
    else{
        return(selected_sets)
    }
}

getXglmnet <- function(x, clusters, n_clusters, type, prototypes=NA){
    # Creates design matrix for glmnet by dealing with clusters (for
    # type="protolasso", discards all cluster members except prototype; for
    # type="clusterRepLasso", replaces all cluster members with a simple
    # average of all the cluster members).

    # Check inputs

    stopifnot(is.matrix(x))
    n <- nrow(x)
    p <- ncol(x)

    stopifnot(is.list(clusters))
    stopifnot(all(lengths(clusters) >= 1))

    stopifnot(length(n_clusters) == 1)
    stopifnot(is.integer(n_clusters) | is.numeric(n_clusters))
    stopifnot(n_clusters >= 0)
    if(type=="protolasso"){
        stopifnot(n_clusters == length(prototypes))
    }
    stopifnot(n_clusters == length(clusters))

    stopifnot(length(type) == 1)
    stopifnot(is.character(type))
    stopifnot(type %in% c("protolasso", "clusterRepLasso"))

    if(type=="protolasso"){
        stopifnot(!is.na(prototypes))
        stopifnot(is.integer(prototypes))
        # if(is.list(clusters)){
            # stopifnot(length(prototypes) == length(clusters))
        # } else{
        #     stopifnot(length(prototypes) == 1)
        # }
        stopifnot(all(!is.na(prototypes)))
        stopifnot(length(prototypes) == length(unique(prototypes)))
    }

    if(n_clusters > 0){
        # if(is.list(clusters)){
        for(i in 1:n_clusters){
            # X_cluster_i <- rowMeans(x[, clusters[[i]]])
            if(type == "protolasso"){
                X_cluster_i <- x[, prototypes[i]]
            } else if(type=="clusterRepLasso"){
                if(length(clusters[[i]]) > 1){
                    X_cluster_i <- rowMeans(x[, clusters[[i]]])
                } else{
                    stopifnot(length(clusters[[i]]) == 1)
                    X_cluster_i <- x[, clusters[[i]]]
                }
                
            } else{
                stop("!(type %in% c(protolasso, clusterRepLasso))")
            }
            
            if(i == 1){
                X_cluster <- as.matrix(X_cluster_i)
                cluster_members <- clusters[[i]]
            } else{
                X_cluster <- cbind(X_cluster, X_cluster_i)
                cluster_members <- c(cluster_members, clusters[[i]])
            }
        }
        non_cluster_feats <- setdiff(1:p, cluster_members)
        # } else{
        #     cluster_members <- clusters
        #     non_cluster_feats <- setdiff(1:p, cluster_members)
        #     if(type == "protolasso"){
        #         X_cluster_i <- x[, prototypes]
        #         stopifnot(all(!(prototypes %in% non_cluster_feats)))
        #     } else if(type=="clusterRepLasso"){
        #         X_cluster_i <- rowMeans(x[, clusters])
        #     } else{
        #         stop("!(type %in% c(protolasso, clusterRepLasso))")
        #     }
        #     X_cluster <- as.matrix(X_cluster_i)
            
        # }
        if(length(non_cluster_feats) > 0){
            X_glmnet <- cbind(X_cluster, x[, non_cluster_feats])
        } else{
            X_glmnet <- X_cluster
        }
        stopifnot(ncol(X_glmnet) == length(non_cluster_feats) + n_clusters)
    } else{
        X_glmnet <- x
        if(type == "protolasso"){stopifnot(length(prototypes) == 0)}
        cluster_members <- integer()
        non_cluster_feats <- 1:p
    }

    colnames(X_glmnet) <- character()

    # Check output
    stopifnot(is.matrix(X_glmnet))
    stopifnot(nrow(X_glmnet) == n)
    stopifnot(ncol(X_glmnet) <= p)
    stopifnot(ncol(X_glmnet) >= 1)

    stopifnot(is.numeric(cluster_members) | is.integer(cluster_members))
    stopifnot(length(cluster_members) <= p)
    stopifnot(length(cluster_members) >= 0)
    if(length(cluster_members) >= 1){
        stopifnot(length(cluster_members) == length(unique(cluster_members)))
        stopifnot(all(round(cluster_members) == cluster_members))
        stopifnot(all(cluster_members) %in% 1:p)
    }

    stopifnot(is.numeric(non_cluster_feats) | is.integer(non_cluster_feats))
    stopifnot(length(non_cluster_feats) <= p)
    stopifnot(length(non_cluster_feats) >= 0)
    if(length(non_cluster_feats) >= 1){
        stopifnot(length(non_cluster_feats) == length(unique(non_cluster_feats)))
        stopifnot(all(round(non_cluster_feats) == non_cluster_feats))
        stopifnot(all(non_cluster_feats) %in% 1:p)
    }
    
    stopifnot(length(intersect(cluster_members, non_cluster_feats)) == 0)
    stopifnot(length(cluster_members) + length(non_cluster_feats) == p)

    return(list(X_glmnet=X_glmnet, cluster_members=cluster_members,
        non_cluster_feats=non_cluster_feats))

}

# protolasso <- function(x, y, R, var_names=NA, nlambda=4000){

#     # R
#     # Numeric p x p matrix; entry ij contains the "substitutive value" of 
#     # feature i for feature j (diagonal must consist of ones, all entries must be between 0 and 1, and matrix must
#     # be symmetric)

#     if(is.data.frame(x)){
#         x <- model.matrix(~ ., x)
#     }

#     stopifnot(is.matrix(x))
#     n <- nrow(x)

#     stopifnot(is.numeric(y) | is.integer(y))
#     stopifnot(n == length(y))

#     p <- ncol(x)

#     stopifnot(is.numeric(y))

#     stopifnot(is.matrix(R))
#     stopifnot(all(dim(R) == p))
#     stopifnot(all(diag(R) == 1))
#     stopifnot(identical(R, t(R)))
#     stopifnot(all(!is.na(R)))
#     stopifnot(all(R >= 0))
#     stopifnot(all(R <= 1))

#     colnames(x) <- character()

#     var_names_provided <- FALSE
#     # If var_names is provided, convert avg_feats entries to character vectors
#     if(all(!is.na(var_names))){
#         stopifnot(all(!is.na(var_names)))
#         var_names_provided <- TRUE
#         stopifnot(length(var_names) == ncol(x))
#     }

#     # Identify clustered features: row i in R contains all the features
#     # that are in a cluster with feature i, so come up with a list containing
#     # all the clusters and then remove the repeats
#     cluster_results <- formatClusters(clusters=NA, R=R, get_prototypes=TRUE,
#         x=x, y=y, p=p)

#     clusters <- cluster_results$clusters
#     # multiple <- cluster_results$multiple
#     prototypes <- cluster_results$prototypes

#     rm(cluster_results)

#     if(is.list(clusters)){
#         n_clusters <- length(clusters)
#     } else{
#         n_clusters <- 1
#     }
    
#     stopifnot(n_clusters == length(prototypes))

#     getXglmnet_results <- getXglmnet(x, clusters, n_clusters,
#         type="protolasso", prototypes=prototypes)

#     X_glmnet <- getXglmnet_results$X_glmnet
#     cluster_members <- getXglmnet_results$cluster_members
#     non_cluster_feats <- getXglmnet_results$non_cluster_feats

#     rm(getXglmnet_results)

#     stopifnot(nrow(X_glmnet) == n)

#     fit <- glmnet(x=X_glmnet, y=y, family="gaussian", alpha=1, nlambda=nlambda)
#     lasso_sets <- unique(predict(fit, type="nonzero"))

#     selected_sets <- getClusterSelsFromGlmnet(lasso_sets, clusters,
#         prototypes, n_cluster_members=length(cluster_members),
#         non_cluster_feats=non_cluster_feats, 
#         var_names_provided=var_names_provided, var_names=var_names,
#         averaging=FALSE)

#     return(list(selected_sets=selected_sets,
#         non_cluster_feats=non_cluster_feats, beta=fit$beta))

# }


# clusterRepLasso <- function(x, y, R, var_names=NA, nlambda=4000){

#     # R
#     # Numeric p x p matrix; entry ij contains the "substitutive value" of 
#     # feature i for feature j (diagonal must consist of ones, all entries must be between 0 and 1, and matrix must
#     # be symmetric)

#     if(is.data.frame(x)){
#         x <- model.matrix(~ ., x)
#     }

#     colnames(x) <- character()

#     stopifnot(is.matrix(x))
#     n <- nrow(x)

#     stopifnot(is.numeric(y) | is.integer(y))
#     stopifnot(n == length(y))

#     p <- ncol(x)

#     stopifnot(is.numeric(y))

#     stopifnot(is.matrix(R))
#     stopifnot(all(dim(R) == p))
#     stopifnot(all(diag(R) == 1))
#     stopifnot(identical(R, t(R)))
#     stopifnot(all(!is.na(R)))
#     stopifnot(all(R >= 0))
#     stopifnot(all(R <= 1))

#     var_names_provided <- FALSE
#     # If var_names is provided, convert avg_feats entries to character vectors
#     if(all(!is.na(var_names))){
#         stopifnot(all(!is.na(var_names)))
#         var_names_provided <- TRUE
#         stopifnot(length(var_names) == ncol(x))
#     }

#     # Identify clustered features: row i in R contains all the features
#     # that are in a cluster with feature i, so come up with a list containing
#     # all the clusters and then remove the repeats

#     # TODO: figure out how to eliminate need for prototypes argument in
#     # getClusterSelsFromGlmnet under clusterRepLasso so that I don't have to get
#     # the prototypes here (shouldn't be necessary for anything)
#     cluster_results <- formatClusters(clusters=NA, R=R, get_prototypes=TRUE,
#         x=x, y=y, p=p)

#     clusters <- cluster_results$clusters
#     # multiple <- cluster_results$multiple
#     prototypes <- cluster_results$prototypes

#     rm(cluster_results)

#     if(is.list(clusters)){
#         n_clusters <- length(clusters)
#     } else{
#         n_clusters <- 1
#     }

#     stopifnot(n_clusters == length(prototypes))

#     getXglmnet_results <- getXglmnet(x, clusters, n_clusters,
#         type="clusterRepLasso")

#     X_glmnet <- getXglmnet_results$X_glmnet
#     cluster_members <- getXglmnet_results$cluster_members
#     non_cluster_feats <- getXglmnet_results$non_cluster_feats

#     rm(getXglmnet_results)

#     stopifnot(nrow(X_glmnet) == n)

#     fit <- glmnet(x=X_glmnet, y=y, family="gaussian", alpha=1, nlambda=nlambda)
#     lasso_sets <- unique(predict(fit, type="nonzero"))

#     cluster_sel_results <- getClusterSelsFromGlmnet(lasso_sets, clusters,
#         prototypes, n_cluster_members=length(cluster_members),
#         non_cluster_feats=non_cluster_feats,
#         var_names_provided=var_names_provided, var_names=var_names,
#         averaging=TRUE)

#     return(list(selected_sets=cluster_sel_results$selected_sets,
#         to_avg_list=cluster_sel_results$to_avg_list,
#         selected_clusts_list=cluster_sel_results$selected_clusts_list,
#         weights_list=cluster_sel_results$weights_list, 
#         non_cluster_feats=non_cluster_feats, beta=fit$beta))
# }

power_law_model <- function(n, p, s, snr, cor, nfolds) {
    # Check inputs
    stopifnot(length(n) == 1)
    stopifnot(n == round(n))
    stopifnot(n >= N_MIN)

    stopifnot(length(s) == 1)
    stopifnot(s == round(s))
    stopifnot(s >= 0)

    stopifnot(length(p) == 1)
    stopifnot(p == round(p))
    stopifnot(p >= s)

    stopifnot(length(snr) == 1)
    stopifnot(is.numeric(snr) | is.integer(snr))
    stopifnot(snr > 0)

    stopifnot(length(cor) == 1)
    stopifnot(is.numeric(cor) | is.integer(cor))
    stopifnot(cor >= 0)
    stopifnot(cor <= 1)

    coefs <- rep(0, p)
    coefs[1:s] <- 1/sqrt(1:s)
    
    my_model <- new_model(name = "power_law_model",
        label = sprintf("Power law model (n= %s, p= %s, s= %s, snr= %s)",
            n, p, s, snr),
        params = list(n = n, p = p, s = s, snr = snr, cor = cor,
            coefs = coefs, nfolds = nfolds),
        simulate = sim_func_selstab
    )

    return(my_model)
}

constant_model <- function(n, p, s, snr, cor, nfolds) {
    # Check inputs
    stopifnot(length(n) == 1)
    stopifnot(n == round(n))
    stopifnot(n >= N_MIN)

    stopifnot(length(s) == 1)
    stopifnot(s == round(s))
    stopifnot(s >= 0)

    stopifnot(length(p) == 1)
    stopifnot(p == round(p))
    stopifnot(p >= s)

    stopifnot(length(snr) == 1)
    stopifnot(is.numeric(snr) | is.integer(snr))
    stopifnot(snr > 0)

    stopifnot(length(cor) == 1)
    stopifnot(is.numeric(cor) | is.integer(cor))
    stopifnot(cor >= 0)
    stopifnot(cor <= 1)

    coefs <- rep(0, p)
    coefs[1:s] <- 1
    
    my_model <- new_model(name = "constant_model",
        label = sprintf("Constant model (n= %s, p= %s, s= %s, snr= %s)",
            n, p, s, snr),
        params = list(n = n, p = p, s = s, snr = snr, cor = cor,
            coefs = coefs, nfolds = nfolds),
        simulate = sim_func_selstab
    )

    return(my_model)
}

sim_func_selstab <- function(n, p, s, snr, cor, coefs, nsim){
    # List we'll return: will be length nsim, and every element will be a
    # named list with elements X and y.
    ret_list <- list()
    for(i in 1:nsim){

        # Generate one draw of mu, X, z, sd, y
        gen_x_mu_res <- gen_x_mu(n, p, s, snr, cor, coefs)
        
        y <- gen_x_mu_res$mu + gen_x_mu_res$sd * rnorm(n)

        stopifnot(is.numeric(y))
        stopifnot(length(y) == n)

        ret_list[[i]] <- list(X = gen_x_mu_res$X, y = y)
    }

    # Check output
    stopifnot(length(ret_list) == nsim)
    for(i in 1:nsim){
        stopifnot(all(names(ret_list[[i]]) == c("X", "y")))
    }
    
    return(ret_list)
}

gen_x_mu <- function(n, p, s, snr, cor, coefs){
    # Check inputs
    stopifnot(length(n) == 1)
    stopifnot(n == round(n))
    stopifnot(n >= N_MIN)

    stopifnot(length(s) == 1)
    stopifnot(s == round(s))
    stopifnot(s >= 0)

    stopifnot(length(p) == 1)
    stopifnot(p == round(p))
    stopifnot(p >= s)

    stopifnot(length(snr) == 1)
    stopifnot(snr > 0)

    stopifnot(length(cor) == 1)
    stopifnot(cor >= 0)
    stopifnot(cor <= 1)

    stopifnot(length(coefs) == p)
    stopifnot(sum(coefs != 0) == s)
    stopifnot(sum(coefs == 0) == p - s)

    # Covariance matrix
    Sigma <- diag(p)
    for(i in 1:round(s/2)){
        for(j in 1:round(s/2)){
            Sigma[i, j] <- Sigma[i, j] + as.numeric(i != j)*cor
        }
    }

    if(cor != 0){
        # Determine number of the p - s noise features that will be in a cluster
        # with each of the s signal features
        cluster_size <- floor((p - s)/s)

        for(i in 1:s){
            for(j in (s + cluster_size*(i - 1) + 1):(s + cluster_size*i)){
                Sigma[i, j] <- Sigma[i, j] + cor
                Sigma[j, i] <- Sigma[j, i] + cor
            }

            for(j in (s + cluster_size*(i - 1) + 1):(s + cluster_size*i)){
                for(k in (s + cluster_size*(i - 1) + 1):(s + cluster_size*i)){
                    Sigma[j, k] <- Sigma[j, k] + as.numeric(j != k)*cor
                }
            }
        }

    }

    x <- mvrnorm(n = n, mu = rep(0, p), Sigma = Sigma)

    mu <- as.numeric(x %*% coefs)

    sd <- sqrt(sum(mu^2) / (n * snr)) # taking snr = (||mu||^2/n) / sigma^2

    # Check output

    stopifnot(is.numeric(x))
    stopifnot(all(dim(x) == c(n, p)))

    stopifnot(is.numeric(mu))
    stopifnot(length(mu) == n)

    stopifnot(length(sd) == 1)
    stopifnot(is.numeric(sd))
    stopifnot(sd > 0)

    return(list(X=x, mu=mu, sd=sd))
}

lasso_cv_selstab <- new_method("lasso_cv_selstab", "Lasso (cross-validation)",
    method = function(model, draw) {

    # Check inputs
    stopifnot(length(names(draw)) == 2)
    stopifnot(all(names(draw) == c("X", "y")))

    stopifnot(is.matrix(draw$X))
    stopifnot(is.numeric(draw$X))

    p <- ncol(draw$X)
    n <- nrow(draw$X)

    stopifnot(is.numeric(draw$y))
    stopifnot(length(draw$y) == n)

    # Use LOOCV (nfolds = n)
    fit <- cv.glmnet(x=draw$X, y=draw$y, nfolds=model$nfolds, grouped=FALSE)

    # Choose selected set by cross-validated MSE
    selected <- predict(fit, type="nonzero", s="lambda.min")$lambda.min

    # Check output
    stopifnot(is.integer(selected))
    stopifnot(length(selected) <= p)
    stopifnot(all(selected %in% 1:p))

    return(list("selected" = selected))

    }
)

get_stable_lambda <- function(X, y, lambdas, nfolds){

    p <- ncol(X)
    n <- nrow(X)

    n_lambdas <- length(lambdas)

    stopifnot(n_lambdas >= 1)

    # Data.frame that will contain selected matrices

    sel_mats <- list()

    for(j in 1:n_lambdas){
        sel_mats[[j]] <- matrix(as.integer(NA), nfolds, p)
    }

    stopifnot(length(sel_mats) == n_lambdas)

    # TODO(gfaletto): change for loop below to allow for k-fold cross-validation
    # (currently hard-coded for LOOCV)
    stopifnot(nfolds == n)
    # Cross-validation
    for(i in 1:nfolds){
        X_i <- X[-i, ]
        stopifnot(all(dim(X_i) == c(n - 1, p)))

        y_i <- y[-i]
        stopifnot(length(y_i) == n - 1)

        fit_i<- glmnet(x=X_i, y=y_i, lambda=lambdas)

        coefs <- as.matrix(coef(fit_i))
        stopifnot(all(dim(coefs) == c(p + 1, n_lambdas)))

        # Eliminate intercept row
        coefs <- coefs[2:(p + 1), ]

        # Convert to indicators of observed features
        coefs <- (coefs != 0)*1

        # Store in selection matrices
        for(j in 1:n_lambdas){
            stopifnot(all(is.na(sel_mats[[j]][i, ])))

            sel_mats[[j]][i, ] <- as.integer(coefs[, j])
        }
    }

    # Check selection matrices
    for(j in 1:n_lambdas){
        stopifnot(all(!is.na(sel_mats[[j]])))
        stopifnot(all(sel_mats[[j]] %in% c(0, 1)))
    }

    # Calculate stability of each matrix of selections

    stabs <- rep(as.numeric(NA), n_lambdas)

    for(j in 1:n_lambdas){
        stopifnot(is.na(stabs[j]))
        stabs[j] <- calcNSBStab(sel_mats[[j]])
    }

    stopifnot(all(!is.na(stabs)))

    if(PRINT_STAB_LAMDDA){
        gg_df <- data.frame("Lambda" = lambdas, "Stability" = stabs)
        plot <- ggplot(gg_df, aes(x=Lambda, y=Stability)) + geom_point()
        print(plot)
    }

    # # Choose largest lambda maximizing stability
    # lambda_star <- max(lambdas[stabs == max(stabs)])

    # Set of lambdas that maximize stability
    lambda_cands <- lambdas[stabs == max(stabs)]

    if(length(lambda_cands) == 1){
        lambda_star <- lambda_cands
    } else{
        # If more than one lambda maximizes stability, choose the one with
        # the lowest cross-validation error
        # lambda_star <- cv.glmnet(X, y, lambda=lambda_cands)$lambda.min
        cvg <- cv.glmnet(X, y, lambda=lambdas)
        # cvg$lambda may differ from lambdas slightly due to rounding, but
        # should be close (and we need to check that it's in the same order
        # in order for code to work)
        stopifnot(all(abs(cvg$lambda - lambdas) < 10^(-6)))

        stopifnot(sum(lambdas %in% lambda_cands) == length(lambda_cands))
        cvg_lambda_cands <- lambdas[lambdas %in% lambda_cands]
        lambda_cands_errors <- cvg$cvm[lambdas %in% lambda_cands]
        lambda_star <- cvg_lambda_cands[which.min(lambda_cands_errors)]
        stopifnot(lambda_star %in% lambda_cands)
    }

    stopifnot(is.numeric(lambda_star))
    stopifnot(length(lambda_star) == 1)
    stopifnot(lambda_star %in% lambdas)

    return(lambda_star)
}

sel_stab <- function(X, y, lambdas, nfolds){

    # Check inputs

    stopifnot(is.matrix(X))
    stopifnot(is.numeric(X))

    p <- ncol(X)
    n <- nrow(X)

    stopifnot(is.numeric(y))
    stopifnot(length(y) == n)

    stopifnot(nfolds == round(nfolds))
    stopifnot(nfolds >= 3)
    stopifnot(nfolds <= n)
    # Use LOOCV (nfolds = n)
    stopifnot(nfolds == n)

    stopifnot(is.numeric(lambdas))

    # Get most (empirically) stable lambda
    lambda_star <- get_stable_lambda(X, y, lambdas, nfolds)

    # Fit lasso
    fit <- glmnet(x=X, y=y, lambda=lambdas)

    pred <- predict(fit, type="nonzero", s=lambda_star)

    if(is.data.frame(pred)){
        selected <- pred[, 1]
    } else{
        selected <- integer(0)
    }
    
    # Check output
    stopifnot(is.integer(selected))
    stopifnot(length(selected) <= p)
    stopifnot(all(selected %in% 1:p))

    return(selected)
}

sel_stab_cv <- function(X, y, lambdas, nfolds){

    # Check inputs

    stopifnot(is.matrix(X))
    stopifnot(is.numeric(X))

    p <- ncol(X)
    n <- nrow(X)

    stopifnot(is.numeric(y))
    stopifnot(length(y) == n)

    stopifnot(nfolds == round(nfolds))
    stopifnot(nfolds >= 3)
    stopifnot(nfolds <= n)

    # For now, always use LOOCV (nfolds = n)
    stopifnot(nfolds == n)

    stopifnot(is.numeric(lambdas))

    # Get optimal lambda for CV MSE
    # TODO: FIX THIS--THIS DOESN'T PROVIDE THE MINIMUM MSE LAMBDA
    fit <- cv.glmnet(x=X, y=y, nfolds=nfolds, grouped=FALSE)
    mse_lambda <- fit$lambda.min

    # Get most (empirically) stable lambda
    lambda_star <- get_stable_lambda(X, y, lambdas, nfolds)

    # Average of these
    mean_lambda <- exp((log(mse_lambda) + log(lambda_star))/2)

    # Add to lambdas
    lambdas <- sort(c(lambdas, mean_lambda), decreasing=TRUE)

    # Fit lasso
    fit <- glmnet(x=X, y=y, lambda=lambdas)

    pred <- predict(fit, type="nonzero", s=mean_lambda)

    if(is.data.frame(pred)){
        selected <- pred[, 1]
    } else{
        selected <- integer(0)
    }
    
    # Check output
    stopifnot(is.integer(selected))
    stopifnot(length(selected) <= p)
    stopifnot(all(selected %in% 1:p))

    return(selected)
}

lasso_stab <- new_method("lasso_stab", "Lasso (stability)",
    method = function(model, draw) {

        # Check inputs

        stopifnot(length(names(draw)) == 2)
        stopifnot(all(names(draw) == c("X", "y")))

        # Get lambda sequence the same way glmnet does
        lambdas <- cv.glmnet(x=draw$X, y=draw$y)$lambda

        stopifnot(is.numeric(lambdas))

        return(list("selected" = sel_stab(draw$X, draw$y, lambdas,
            model$nfolds)))
    }
)

## TODO: create function that gets lambda sequence like glmnet does. (Need to
# crib from function get_start in glmnetFlex.R.)
# getGlmnetLambdas <- function(nobs, nvars, lambda_max=as.numeric(NA),
#     lambda.min.ratio=ifelse(nobs < nvars, 0.01, 1e-04), nlambda=100){
#     # Check inputs

#     stopifnot(length(nlambda) == 1)
#     stopifnot(is.numeric(nlambda) | is.integer(nlmabda))
#     stopifnot(round(nlambda) == nlambda)
#     stopifnot(nlambda >= 1)

#     # compute lambda max
#     ju <- rep(1, nvars)
#     r <- y - mu
#     eta <- family$linkfun(mu)
#     v <- family$variance(mu)
#     m.e <- family$mu.eta(eta)
#     weights <- weights / sum(weights)
#     rv <- r / v * m.e * weights
#     if (inherits(x, "sparseMatrix")) {
#         xm <- attr(x, "xm")
#         xs <- attr(x, "xs")
#         g <- abs((drop(t(rv) %*% x) - sum(rv) * xm) / xs)
#     } else {
#         g <- abs(drop(t(rv) %*% x))
#     }
#     g <- g * ju / ifelse(vp > 0, vp, 1)
#     lambda_max <- max(g) / max(alpha, 1e-3)

#     return(exp(seq(log(lambda_max), log(lambda_max * lambda.min.ratio),
#                 length.out = nlambda)))
# }

lasso_stab_cv <- new_method("lasso_stab_cv", "Lasso (stability + MSE)",
    method = function(model, draw) {

        # Check inputs

        stopifnot(length(names(draw)) == 2)
        stopifnot(all(names(draw) == c("X", "y")))

        lambdas <- cv.glmnet(x=draw$X, y=draw$y)$lambda

        stopifnot(is.numeric(lambdas))

        return(list("selected" = sel_stab_cv(draw$X, draw$y, lambdas,
            model$nfolds)))
    }
)

prop_discs <- new_metric("prop_discs", "Proportion of Discoeveries",
    metric = function(model, out) {

    # Check inputs
    stopifnot(is.numeric(model$coefs))
    stopifnot(length(model$coefs) == model$p)

    stopifnot(is.integer(out$selected))
    stopifnot(length(out$selected) >= 0)
    stopifnot(length(out$selected) <= model$p)

    true_feats <- as.integer(which(model$coefs != 0))

    stopifnot(length(true_feats) == model$s)

    discoveries <- intersect(true_feats, out$selected)

    stopifnot(is.integer(discoveries))
    stopifnot(length(discoveries) <= min(model$s, length(out$selected)))

    return(length(discoveries)/model$s)
  
    }
)

prop_rejects <- new_metric("prop_rejects", "Proportion of Rejections",
    metric = function(model, out) {

    # Check inputs
    stopifnot(is.numeric(model$coefs))
    stopifnot(length(model$coefs) == model$p)

    stopifnot(is.integer(out$selected))
    stopifnot(length(out$selected) >= 0)
    stopifnot(length(out$selected) <= model$p)

    true_feats <- as.integer(which(model$coefs != 0))

    stopifnot(length(true_feats) == model$s)

    rejects <- setdiff(setdiff(1:p, out$selected), true_feats)

    stopifnot(is.integer(rejects))
    stopifnot(length(rejects) <= min(model$p - model$s,
        model$p -length(out$selected)))

    return(length(rejects)/(model$p - model$s))
  
    }
)

get_ranks <- function(props, coefs){
    # Calculates true ranks of features from a model (given the coefficients
    # of the model) as well as the ranks implied by selection proportions for
    # each feature (yielded by e.g. stability selection).

    # Inputs:

    # props
    # A numeric vector of selection proportions (length must be equal to the
    # number of features).

    # coefs
    # A vector of coefficients of features in the true model (length must be
    # equal to the number of features.)

    # Output: A named list containing two vectors

    # est_ranks: The ranks implied by the selection proportions in props (that
    # is, the ranks estimated by the procedure that yielded the proportions).

    # true_ranks: The true ranks implied by the coefficients in coefs.

    # Check inputs

    p <- length(coefs)

    stopifnot(is.numeric(coefs) | is.integer(coefs))

    # Take absolute value of coefficients, because only size of coefficient
    # is important for ranking, not sign
    coefs <- abs(coefs)

    stopifnot(is.numeric(props))
    stopifnot(length(props) == p)
    stopifnot(all(props >= 0) & all(props <= 1))

    # Get ranking of features using selection proportions (highest proportions
    # are highest ranks). (Doing it this way instead of using rank function
    # because I want to preserve ties.)

    est_ranks <- as.integer(factor(-1*props))

    # Get rankings of features from true model
    true_ranks <- rep(as.integer(NA), p)

    # Need to sort coefficients so that proper ranking is assigned later
    unique_coefs <- sort(unique(coefs), decreasing=TRUE)

    n_unique_coefs <- length(unique_coefs)

    stopifnot(n_unique_coefs <= p)
    stopifnot(n_unique_coefs >= 1)

    for(i in 1:n_unique_coefs){
        inds_i <- coefs == unique_coefs[i]
        stopifnot(sum(inds_i) >= 1)
        stopifnot(all(is.na(true_ranks[inds_i])))
        true_ranks[inds_i] <- i
    }

    # Check output

    stopifnot(all(!is.na(est_ranks)))
    stopifnot(length(est_ranks) == p)
    stopifnot(length(unique(est_ranks)) <= p)
    stopifnot(all(est_ranks %in% 1:p))

    stopifnot(all(!is.na(true_ranks)))
    stopifnot(length(true_ranks) == p)
    stopifnot(length(unique(true_ranks)) <= p)
    stopifnot(all(true_ranks %in% 1:p))

    return(list(est_ranks=est_ranks, true_ranks=true_ranks))
}

rank_cor_kendall <- new_metric("rank_cor_kendall",
    "Ranking Accuracy (Kendall)", metric = function(model, out) {

    # Check inputs (get_ranks function will check the validity of
    # out$proportions and model$coefs)

    stopifnot(length(model$p) == 1)
    stopifnot(is.numeric(model$p) | is.integer(model$p))
    stopifnot(round(model$p) == model$p)

    stopifnot(length(model$coefs) == model$p)

    rank_results <- get_ranks(props=out$proportions, coefs=model$coefs)

    # Similarity metric (Kendall's Tau)
    similarity <- Kendall(rank_results$true_ranks, rank_results$est_ranks)$tau

    # Check output
    stopifnot(is.numeric(similarity))
    stopifnot(similarity >= -1 & similarity <= 1)

    return(similarity)

    }
)

rank_cor_spearman <- new_metric("rank_cor_spearman",
    "Ranking Accuracy (Spearman)", metric = function(model, out) {

    # Check inputs (get_ranks function will check the validity of
    # out$proportions and model$coefs)

    stopifnot(length(model$p) == 1)
    stopifnot(is.numeric(model$p) | is.integer(model$p))
    stopifnot(round(model$p) == model$p)

    stopifnot(length(model$coefs) == model$p)

    rank_results <- get_ranks(props=out$proportions, coefs=model$coefs)

    # Similarity metric (Spearman's Rho) (suppress warnings because otherwise
    # if there are ties it will give a warning about the p-value, which we
    # don't care about, being off)
    similarity <- suppressWarnings(cor.test(rank_results$true_ranks,
        rank_results$est_ranks, method="spearman")$estimate)

    # Check output
    stopifnot(is.numeric(similarity))
    stopifnot(similarity >= -1 & similarity <= 1)

    return(similarity)

    }
)

stabsel_props <- new_method("stabsel_props", "Stability Selection",
    method = function(model, draw) {

        # Check inputs

        stopifnot(length(names(draw)) == 2)
        stopifnot(all(names(draw) == c("X", "y")))
        stopifnot(ncol(draw$X) == model$p)
        stopifnot(nrow(draw$X) == model$n)

        B <- 100
        R <- diag(model$p)

        stopifnot(length(draw$y) == model$n)

        # Determine lambda: do lasso with cross-validation on full data sample.
        # Cross-validated model size will be the size we use for stability
        # selection.
        inds_size <- sample(1:model$n, floor(model$n/2))

        size_results <- cv.glmnet(x=draw$X[inds_size, ], y=draw$y[inds_size],
            parallel=TRUE, family="gaussian")

        fit <- css(X=draw$X, y=draw$y, lambda=size_results$lambda.min,
            R=R, B=B, sampling_type="SS")

        props <- fit$feat_sel_props

        # Check output

        stopifnot(is.numeric(props))
        stopifnot(length(props) == model$p)
        stopifnot(all(props >= 0 & props <= 1))

        return(list(proportions = props))
    }
)

# stabsel_props_add <- new_method("stabsel_props_add",
#     "Stability Selection (add feats)",
#     method = function(model, draw, prop_model_add) {

#         # Check inputs

#         stopifnot(length(names(draw)) == 2)
#         stopifnot(all(names(draw) == c("X", "y")))
#         stopifnot(ncol(draw$X) == model$p)
#         stopifnot(nrow(draw$X) == model$n)

#         B <- 100
#         R <- diag(model$p)

#         stopifnot(length(draw$y) == model$n)

#         stopifnot(length(prop_model_add) == 1)
#         stopifnot(is.numeric(prop_model_add))
#         stopifnot(prop_model_add >= 0)

#         # Determine lambda: do lasso with cross-validation on full data sample.
#         # Cross-validated model size will be the size we use for stability
#         # selection.
#         inds_size <- sample(1:model$n, floor(model$n/2))

#         size_results <- cv.glmnet(x=draw$X[inds_size, ], y=draw$y[inds_size],
#             parallel=TRUE, family="gaussian")

#         # Determine proportion of features to add on every subsample: take
#         # prop_model_add of the features in the true model
#         sel_prop <- prop_model_add*model$s/model$p

#         fit <- css(x=draw$X, y=draw$y, lambda=size_results$lambda.min,
#             R=R, B=B, sampling_type="SS", prop_feats_add=sel_prop)

#         props <- fit$feat_sel_props

#         # Check output

#         stopifnot(is.numeric(props))
#         stopifnot(length(props) == model$p)
#         stopifnot(all(props >= 0 & props <= 1))

#         return(list(proportions = props))
#     },
#     settings = list(prop_model_add=0.25)
# )




getBinMat <- function(output, meth, model_size){
    stopifnot(length(output) == n_sims)

    ret <- matrix(0, n_sims, p)

    for(j in 1:n_sims){
        output_j <- output[[j]]
        cssr_meth <- FALSE
        # Only need models of size model_size from method meth
        if(meth %in% c("elastic_net", "lasso_random")){
            feat_list <- output_j$lasso_selected
        } else if(meth %in% c("clusRepLasso_cssr", "clusRepLasso_cssr_est")){
            feat_list <- output_j$selected_clusts_list
        } else if(meth %in% c("protolasso_cssr", "protolasso_cssr_est")){
            feat_list <- output_j$selected_sets
        } else if(meth %in% c("SS_CSS_avg_cssr", "SS_CSS_sparse_cssr",
            "SS_CSS_weighted_cssr", "SS_CSS_avg_cssr_est",
            "SS_CSS_sparse_cssr_est", "SS_CSS_weighted_cssr_est",
            "SS_SS_cssr")){
            feat_list <- output_j$selected
            cssr_meth <- TRUE
        }

        if(!cssr_meth){
            feat_ind <- lengths(feat_list) == model_size
            if(any(feat_ind)){
                feats_j <- feat_list[[min(which(feat_ind))]]
                if(is.list(feats_j)){
                    feats_j <- unlist(feats_j)
                }
                ret[j, feats_j] <- 1
            }
        } else{
            if(length(feat_list) >= model_size){
                feats_j <- feat_list[[model_size]]
                if(!is.null(feats_j)){
                    ret[j, feats_j] <- 1
                }
            }
            
        }
    }
    return(ret)
}



createLossesPlot3 <- function(df_gg, n_methods, legend=TRUE,
    plot_errors=TRUE, subtitle=FALSE){

    if((sig_blocks + k_unblocked) %% 2 == 0){
        max_rank <- sig_blocks + k_unblocked
    } else{
        max_rank <- sig_blocks + k_unblocked + 1
    }

    plot <- ggplot(df_gg, aes(x=ModelSize, y=MSE, color=Method,
        shape=Method)) + scale_shape_manual(values=1:n_methods) +
        suppressWarnings(geom_point(size=2.5, alpha=1)) + 
        xlab("No. Fitted Coefficients") +
        scale_x_continuous(breaks=seq(2, max_rank, by=2))

    if(subtitle){
        subtitle_txt <- paste("n = ", n_model, ", p = ", p, ",
            beta_high = ", beta_high, sep="")
        plot <- plot + labs(subtitle=subtitle_txt)
    }

    if(plot_errors){
        plot <- plot + geom_errorbar(aes(ymin = MSELower, ymax = MSEUpper),
            width = 0.5)
    }

    if(!legend){
        plot <- plot + theme(legend.position="none")
    }

    return(plot)
}

createNSBStabPlot2 <- function(df_gg, legend=TRUE, plot_errors=TRUE, 
    subtitle=FALSE){
    require(ggplot2)
    
    # Get max ranking value for integer labels on horizontal axis
    max_rank <- max(df_gg$ModelSize)
    if((max_rank %% 2) != 0){
        max_rank <- max_rank + 1
    }

    if(subtitle){
        subtitle_txt <- paste("n = ", n_model, ", p = ", p, ",
            beta_high = ", beta_high, sep="")
    }

    plot <- ggplot(df_gg, aes(x=ModelSize, y=NSBStability,
        color=Method, shape=Method)) + scale_shape_manual(values=1:n_methods)

    plot <- plot + suppressWarnings(geom_point(size=2.5, alpha=1)) +
        xlab("No. Fitted Coefficients") + ylab("NSB Stability") +
        scale_x_continuous(breaks=seq(2, max_rank, by=2))

    if(subtitle){
        plot <- plot + labs(subtitle=subtitle_txt)
    }

    if(plot_errors){
        plot <- plot + geom_errorbar(aes(ymin = StabLower, ymax = StabUpper),
            width = 0.5)
    }

    if(!legend){
        plot <- plot + theme(legend.position="none")
    }

    return(plot)
}




createStabMSEPlot2 <- function(df_gg, n_methods, legend=TRUE, plot_errors=FALSE,
    subtitle=FALSE){

    require(ggplot2)

    if(subtitle){
        subtitle_txt <- paste("n = ", n_model, ", p = ", p, ", k = ", k, ",
            beta_high = ", beta_high, sep="")
    }

    plot <- ggplot(df_gg, aes(x=NSBStability, y=MSE,
        color=Method, shape=Method)) + scale_shape_manual(values=1:n_methods) +
        suppressWarnings(geom_point(size=2.5, alpha=1)) + xlab("NSB Stability")
        # + 
        # geom_hline(yintercept=true_model_loss
        #     # , color="red"
        #     , linetype = "dashed"
        #     )

    if(subtitle){
        plot <- plot + labs(subtitle=subtitle_txt)
    }

    if(plot_errors){
        plot <- plot + geom_errorbar(aes(ymin = MSELower, ymax = MSEUpper)
            , width = 0.05
            )
        if(cluster_count == "none"){
            plot <- plot + geom_errorbar(aes(xmin = StabLower, xmax = StabUpper)
            , width = 0.25
            )
        }
    }

    if(!legend){
        plot <- plot + theme(legend.position="none")
    }

    return(plot)
}


createPhatPlot2 <- function(output, ylab="Proportion of Subsamples",
    feats_to_display=20, line=FALSE, tau=NA, title=FALSE){
    # Creates bar plot of estimated probability of selection under our
    # proposed modification to stability selection against ranking of these
    # probabilities, in descending order of ranking (high to low estimated
    # probabilities). Bars are color coded according to whether feature
    # is blocked or unblocked, signal or noise.

    # Outputs: ggplot2 object of plot

    require(ggplot2)

    if(line & is.na(tau)){
        stop("If line is TRUE, must provide tau")
    }

    # p <- sim_model@params$p

    # beta_high <- sim_model@params$beta_high

    stopifnot(length(output) == n_sims)

    output_1 <- output[[1]]

    sel_mat <- output_1$css_res$feat_sel_mat

    stopifnot(p == ncol(sel_mat))

    props <- colMeans(sel_mat)

    stopifnot(nblocks == 1)
    stopifnot(sig_blocks == 1)

    feat_names <- rep("Noise Features", p)
    # Overwrite with other feature names as appropriate
    feat_names[1:block_size] <- "Proxies For Z"
    feat_names[(block_size + 1):(block_size + k_unblocked)] <- "Weak Signal Features"

    stopifnot(length(props) == p)
    stopifnot(length(feat_names) == p)
    stopifnot(length(1:p) == p)

    data <- data.frame(Feature=1:p, props, phat=props, Legend=feat_names)

    data <- drop(data)

    # Sort data in descending order of p_hat
    data <- data[order(data$phat, decreasing=TRUE), ]

    # Add a ranking variable
    colnames_orig <- colnames(data)
    data <- data.frame(1:p, data)
    colnames(data) <- c("Ranking", colnames_orig)

    # Keep only first feats_to_display features (easier to visualize)
    data <- data[1:min(nrow(data), feats_to_display), ]

    subtitle <- paste("p = ", p ,", beta_high = ", beta_high, sep="") 

    plot <- ggplot(data, aes(x=Ranking, y=phat, fill=Legend)) + geom_col() +
        ylab(ylab) 

    if(title){
        plot <- plot + labs(subtitle=subtitle)
    }

    if(line){
        plot <- plot + geom_hline(yintercept=tau, linetype="dashed",
            color="black") 
    }

    print(plot)

    return(plot)
}

saveFigure2 <- function(subdir, plot, size="large", filename){

    stopifnot(length(filename) == 1)

    stopifnot(size %in% c("small", "medium", "mlarge", "large", "xlarge",
        "slide", "xmlarge"))

    dir_fig <- file.path(wd, subdir)
    dir.create(dir_fig, showWarnings = FALSE, recursive = TRUE)
    setwd(dir_fig)

    if(size=="large"){
        pdf(file=filename, title=filename, paper="special", width=10.5,
            height=5.5, bg="white")
    } else if(size=="mlarge"){
        pdf(file=filename, title=filename, paper="special", width=10.5,
            height=4, bg="white")
    } else if(size=="xlarge"){
        pdf(file=filename, title=filename, paper="special", width=10.5,
            height=6.5, bg="white")
    } else if(size=="xmlarge"){
        pdf(file=filename, title=filename, paper="special", width=9,
            height=6.5, bg="white")
    } else if(size=="small"){
        pdf(file=filename, title=filename, paper="special", width=5, height=5,
            bg="white")
    } else if(size=="medium"){
        pdf(file=filename, title=filename, paper="special", width=6.4,
            height=5.5, bg="white")
    } else if(size=="slide"){
        pdf(file=filename, title=filename, paper="special", width=540/72,
            height=360/72, bg="white")
    }
    
    plot(plot)
    dev.off()

    setwd(wd)
}

genPlotDf <- function(completed_sim, alpha=0.05){
    print("getting evals...")
    e <- evals(completed_sim)
    print("done! converting to data.frame...")
    t0 <- Sys.time()
    edf <- as.data.frame(e)
    print("done! time to convert to data.frame:")
    print(Sys.time() - t0)
    stopifnot("cssr_mse" %in% colnames(edf))
    print("Getting MSEs and stability metrics...")
    t1 <- Sys.time()

    methods <- unique(edf$Method)
    n_methods <- length(methods)

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
        stopifnot(n_sims > 1)

        meth_i_vec <- rep(as.numeric(NA), sig_blocks + k_unblocked)
        o_i <- output(completed_sim, methods=methods[i])@out

        for(k in 1:(sig_blocks + k_unblocked)){
            # MSE
            inds_k <- (sig_blocks + k_unblocked)*(0:(n_sims - 1)) + k
            stopifnot(all(inds_k %in% 1:nrow(edf_i)))
            stopifnot("cssr_mse" %in% colnames(edf_i))
            mses_ik <- edf_i[inds_k, "cssr_mse"]

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
        print(paste("Done with method", methods[i]))
        print("Time in this step so far:")
        t_i_end <- Sys.time()
        print(t_i_end - t_1)
        print("Estimated time to go:")
        print((n_methods - i)/i*(t_i_end - t_1))
    }

    print("Done! time to get metrics:")
    print(Sys.time() - t1)

    results_df <- data.frame(ModelSize=model_sizes, Method=methods_vec, MSE=mses,
        MSELower=mses - margins, MSEUpper=mses + margins, NSBStability=nsbstabs,
        StabLower=nsb_lowers, StabUpper=nsb_uppers)

    return(list(results_df=results_df, n_methods=n_methods))
}