# setwd("/Users/gregfaletto/Google Drive/Data Science/LaTeX/Generalized Stability Selection Presentation")


library(ggplot2)

nameMap <- function(sys_name){
    # Takes in computer-friendly name for method and outputs display name
    if(!is.character(sys_name)){
        stop("!is.character(sys_name)")
    }
    ret <- sys_name
    ret[sys_name %in%  c("lasso", "lasso_random")] <- "Lasso"
    ret[sys_name %in%  c("lassoSS_phat", "SS_SS_random")] <-
        "Stability Selection"
    ret[sys_name %in%  c("lassoSS_phat_ideal", "SS_GSS_random")] <-
        "Generalized Stability Selection"
    ret[sys_name %in%  c("SS_GSS_random_avg")] <-
        "GSS (Weighted Averaging)"
    ret[sys_name %in%  c("SS_GSS_random_avg_unwt")] <-
        "GSS (Unweighted Averaging)"
    ret[sys_name %in%  c("BRVZ_avg_unwt")] <-
        "Lasso on Cluster Average"
    ret[sys_name %in%  c("lasso_proto")] <-
        "Lasso on Cluster Prototype"
    return(ret)
}

makeDfLabelsDisplay <- function(df){
    # Takes in data.frame with computer-friendly labels and returns 
    # data.frame with display-friendly labels
    if(!("Label" %in% colnames(df))){
        stop("!(Label %in% colnames(df))")
    }
    if(ncol(df) != 3){
        stop("ncol(df) != 3 (before modifications)")
    }
    labels <- as.factor(nameMap(as.character(df$Label)))
    labels <- droplevels(labels)
    df$Label <- labels
    if(!all(levels(df$Label) %in% c("Lasso", "Stability Selection",
        "Generalized Stability Selection", "GSS (Weighted Averaging)",
        "GSS (Unweighted Averaging)", "Lasso on Cluster Average",
        "Lasso on Cluster Prototype"))){
        stop("!all(levels(df$Label) %in% c(Lasso, Stability Selection,vGeneralized Stability Selection, GSS (Weighted Averaging), GSS (Unweighted Averaging, Lasso on Cluster Average, Lasso on Cluster Prototype)))")
    }
    if(ncol(df) != 3){
        stop("ncol(df) != 3 (after modifications)")
    }
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

    # print("sim_evals@evals:")
    # print(sim_evals@evals)

    # print("sim_evals@evals[[method]]$r1.1:")
    # print(sim_evals@evals[[method]]$r1.1)

    # print("p:")
    # print(p)
    # print("method:")
    # print(method)
    # print("length(sim_evals@evals[[method]]$r1.1$phat):")
    # print(length(sim_evals@evals[[method]]$r1.1$phat))
    # print("length(sim_evals@evals[[method]]$r1.1$labels):")
    # print(length(sim_evals@evals[[method]]$r1.1$labels))

    if(n_iters > 1){
        # Create list of plots
        plots <- list()
    }

    for(i in 1:n_iters){
        # Index corresponding to this i
        id <- paste("r1.", i, sep="")

        data <- data.frame(1:p, sim_evals@evals[[method]][[id]]$phat,
            sim_evals@evals[[method]][[id]]$labels)

        colnames(data) <- c("Feature", "phat", "Label")

        # re-name noise features

        levels(data$Label) <- c(levels(data$Label), "Noise Features",
            "Proxies For Z", "Weak Signal Features")

        data$Label[data$Label == "Unblocked Noise Features"] <- "Noise Features"
        data$Label[data$Label == "Blocked Noise Features"] <- "Proxies For Z"
        data$Label[data$Label == "Blocked Signal Features"] <- "Proxies For Z"
        data$Label[data$Label == "Unblocked Signal Features"] <- "Weak Signal Features"

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

        plot <- ggplot(data, aes(x=Ranking, y=phat, fill=Label)) + geom_col() +
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

    colnames(data) <- c("Feature", "phat", "Label")

    # re-name noise features

    levels(data$Label) <- c(levels(data$Label), "Noise Features",
        "Proxies For Z", "Weak Signal Features")

    data$Label[data$Label == "Unblocked Noise Features"] <- "Noise Features"
    data$Label[data$Label == "Blocked Noise Features"] <- "Proxies For Z"
    data$Label[data$Label == "Blocked Signal Features"] <- "Proxies For Z"
    data$Label[data$Label == "Unblocked Signal Features"] <- "Weak Signal Features"

    data <- drop(data)

    # Sort data in descending order of p_hat
    data <- data[order(data$phat, decreasing=TRUE), ]

    # Add a ranking variable
    colnames_orig <- colnames(data)
    data <- data.frame(1:p, data)
    colnames(data) <- c("Ranking", colnames_orig)

    # Keep only first feats_to_display features (easier to visualize)
    data <- data[1:min(nrow(data), feats_to_display), ]
 
    plot <- ggplot(data, aes(x=Ranking, y=phat, fill=Label)) + geom_col() +
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
    
    colnames(data) <- c("Feature", "phat", "Label")
    
    # re-name noise features
    
    levels(data$Label) <- c(levels(data$Label), "Noise Features",
                            "Proxies For Z", "Weak Signal Features")
    
    data$Label[data$Label == "Unblocked Noise Features"] <- "Noise Features"
    data$Label[data$Label == "Blocked Noise Features"] <- "Proxies For Z"
    data$Label[data$Label == "Blocked Signal Features"] <- "Proxies For Z"
    data$Label[data$Label == "Unblocked Signal Features"] <- "Weak Signal Features"
    
    data <- drop(data)
    
    # Sort data in descending order of p_hat
    data <- data[order(data$phat, decreasing=TRUE), ]
    
    # Add a ranking variable
    colnames_orig <- colnames(data)
    data <- data.frame(1:p, data)
    colnames(data) <- c("Ranking", colnames_orig)

    # Keep only first feats_to_display features (easier to visualize)
    data <- data[1:min(nrow(data), feats_to_display), ]
    
    plot <- ggplot(data, aes(x=Ranking, y=phat, fill=Label)) + geom_col() +
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

    colnames(data) <- c("Feature", "phat", "Label")

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
    data <- makeDfLabelsDisplay(data)

    # Get max ranking value for integer labels on horizontal axis
    max_rank <- max(data$Ranking)
    if((max_rank %% 2) != 0){
        max_rank <- max_rank + 1
    }

    plot <- ggplot(data) + geom_col(mapping=aes(x=Ranking, y=phat, fill=Label)) +
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
    if(length(J) > 1){
        stop("length(J) > 1")
    }
    if(is.na(J)){
        J <- (nrow(false_selecs) - 1)/n_methods
    }
    if(J != round(J)){
        print(nrow(false_selecs))
        print(n_methods)
        print(J)
        stop("J != round(J)")
    }
    if(any(false_selecs < 0)){
        stop("any(false_selecs < 0)")
    }
    if(n_sims != ncol(j_choices_mat) | nrow(false_selecs) != nrow(j_choices_mat)){
        print(dim(false_selecs))
        print(dim(j_choices_mat))
        stop("ncol(false_selecs) != ncol(j_choices_mat) | nrow(false_selecs) != nrow(j_choices_mat)")
    }

    if(n_sims > 1){
        row_inds <- which(rowSums(j_choices_mat) > 0)
        n_false_selecs <- sum(rowSums(j_choices_mat) > 0)
        false_selecs_vec <- numeric(n_false_selecs)
        j_choices <- which(rowSums(j_choices_mat[1:J, ]) > 0)
        if(length(j_choices)*n_methods != n_false_selecs){
            stop("length(j_choices)*n_methods + 1 != n_false_selecs")
        }
        if(length(row_inds) != n_false_selecs){
            stop("length(row_inds) != n_false_selecs")
        }


        for(i in 1:n_false_selecs){
            cols_i <- j_choices_mat[row_inds[i], ]
            false_selecs_vec[i] <- mean(false_selecs[row_inds[i], cols_i])
        }
    } else{
        row_inds <- which(j_choices_mat)
        n_false_selecs <- sum(j_choices_mat)
        false_selecs_vec <- numeric(n_false_selecs)
        j_choices <- which(j_choices_mat[1:J, 1])
        if(length(j_choices)*n_methods != n_false_selecs){
            stop("length(j_choices)*n_methods + 1 != n_false_selecs")
        }
        if(length(row_inds) != n_false_selecs){
            stop("length(row_inds) != n_false_selecs")
        }

        false_selecs_vec <- false_selecs[row_inds, 1]
    }

    df_gg <- data.frame(rep(j_choices, n_methods), false_selecs_vec,
        rep(names, each=length(j_choices)))

    # # Data.frame for plot
    # df_gg <- data.frame(rep(1:k_ss, n_ss_models), as.vector(false_selecs),
    #     rep(names, each=k))

    colnames(df_gg) <- c("Model_Size", "False_Selections", "Label")

    subtitle <- paste("n = ", n_model, ", p = ", p, ", k = ", k,
        ", beta_high = ", beta_high, sep="")

    # Make labels display-friendly
    df_gg <- makeDfLabelsDisplay(df_gg)

    # Get max ranking value for integer labels on horizontal axis
    max_rank <- max(df_gg$Model_Size)
    if((max_rank %% 2) != 0){
        max_rank <- max_rank + 1
    }

    plot <- ggplot(df_gg, aes(x=Model_Size, y=False_Selections, color=Label,
        shape=Label)) + geom_point(size=2.5, alpha=1) + xlab("Model Size") +
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
    if(length(J) > 1){
        stop("length(J) > 1")
    }
    if(is.na(J)){
        J <- (nrow(proxy_selecs) - 1)/n_methods
    }
    if(J != round(J)){
        print(nrow(proxy_selecs))
        print(n_methods)
        print(J)
        stop("J != round(J)")
    }
    if(any(proxy_selecs < 0)){
        stop("any(proxy_selecs < 0)")
    }
    if(n_sims != ncol(j_choices_mat) | nrow(proxy_selecs) != nrow(j_choices_mat)){
        print(dim(proxy_selecs))
        print(dim(j_choices_mat))
        stop("ncol(proxy_selecs) != ncol(j_choices_mat) | nrow(proxy_selecs) != nrow(j_choices_mat)")
    }

    if(n_sims > 1){
        row_inds <- which(rowSums(j_choices_mat) > 0)
        n_proxy_selecs <- sum(rowSums(j_choices_mat) > 0)
        proxy_selecs_vec <- numeric(n_proxy_selecs)
        j_choices <- which(rowSums(j_choices_mat[1:J, ]) > 0)
        if(length(j_choices)*n_methods != n_proxy_selecs){
            stop("length(j_choices)*n_methods + 1 != n_proxy_selecs")
        }
        if(length(row_inds) != n_proxy_selecs){
            stop("length(row_inds) != n_proxy_selecs")
        }


        for(i in 1:n_proxy_selecs){
            cols_i <- j_choices_mat[row_inds[i], ]
            proxy_selecs_vec[i] <- mean(proxy_selecs[row_inds[i], cols_i])
        }
    } else{
        row_inds <- which(j_choices_mat)
        n_proxy_selecs <- sum(j_choices_mat)
        proxy_selecs_vec <- numeric(n_proxy_selecs)
        j_choices <- which(j_choices_mat[1:J, 1])
        if(length(j_choices)*n_methods != n_proxy_selecs){
            stop("length(j_choices)*n_methods + 1 != n_proxy_selecs")
        }
        if(length(row_inds) != n_proxy_selecs){
            stop("length(row_inds) != n_proxy_selecs")
        }

        proxy_selecs_vec <- proxy_selecs[row_inds, 1]
    }

    df_gg <- data.frame(rep(j_choices, n_methods), proxy_selecs_vec,
        rep(names, each=length(j_choices)))

    # # Data.frame for plot
    # df_gg <- data.frame(rep(1:k_ss, n_ss_models), as.vector(proxy_selecs),
    #     rep(names, each=k))

    colnames(df_gg) <- c("Model_Size", "Proxy_Selections", "Label")

    # Make labels display-friendly
    df_gg <- makeDfLabelsDisplay(df_gg)

    # Get max ranking value for integer labels on horizontal axis
    max_rank <- max(df_gg$Model_Size)
    if((max_rank %% 2) != 0){
        max_rank <- max_rank + 1
    }

    subtitle <- paste("n = ", n_model, ", p = ", p, ", k = ", k, ", beta_high = ",
        beta_high, sep="")

    plot <- ggplot(df_gg, aes(x=Model_Size, y=Proxy_Selections, color=Label,
        shape=Label)) + geom_point(size=2.5, alpha=1) + xlab("Model Size") +
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
    if(length(J) > 1){
        stop("length(J) > 1")
    }
    if(is.na(J)){
        J <- (nrow(proxy_selecs) - 1)/n_methods
    }
    if(J != round(J)){
        print(nrow(proxy_selecs))
        print(n_methods)
        print(J)
        stop("J != round(J)")
    }
    if(any(proxy_selecs < 0)){
        stop("any(proxy_selecs < 0)")
    }
    if(n_sims != ncol(j_choices_mat) | nrow(proxy_selecs) != nrow(j_choices_mat)){
        print(dim(proxy_selecs))
        print(dim(j_choices_mat))
        stop("ncol(proxy_selecs) != ncol(j_choices_mat) | nrow(proxy_selecs) != nrow(j_choices_mat)")
    }

    if(n_sims > 1){
        row_inds <- which(rowSums(j_choices_mat) > 0)
        n_proxy_selecs <- sum(rowSums(j_choices_mat) > 0)
        n_false_selecs <- sum(rowSums(j_choices_mat) > 0)
        proxy_selecs_vec <- numeric(n_proxy_selecs)
        false_selecs_vec <- numeric(n_false_selecs)
        j_choices <- which(rowSums(j_choices_mat[1:J, ]) > 0)
        if(length(j_choices)*n_methods != n_proxy_selecs){
            stop("length(j_choices)*n_methods + 1 != n_proxy_selecs")
        }
        if(length(row_inds) != n_proxy_selecs){
            stop("length(row_inds) != n_proxy_selecs")
        }
        if(length(j_choices)*n_methods != n_false_selecs){
            stop("length(j_choices)*n_methods + 1 != n_false_selecs")
        }
        if(length(row_inds) != n_false_selecs){
            stop("length(row_inds) != n_false_selecs")
        }

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
        if(length(j_choices)*n_methods != n_proxy_selecs){
            stop("length(j_choices)*n_methods + 1 != n_proxy_selecs")
        }
        if(length(row_inds) != n_proxy_selecs){
            stop("length(row_inds) != n_proxy_selecs")
        }
        if(length(j_choices)*n_methods != n_false_selecs){
            stop("length(j_choices)*n_methods + 1 != n_false_selecs")
        }
        if(length(row_inds) != n_false_selecs){
            stop("length(row_inds) != n_false_selecs")
        }

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

    colnames(df_gg) <- c("Model_Size", "Weak_Signal_Selections", "Label")

    # Make labels display-friendly
    df_gg <- makeDfLabelsDisplay(df_gg)

    # Get max ranking value for integer labels on horizontal axis
    max_rank <- max(df_gg$Model_Size)
    if((max_rank %% 2) != 0){
        max_rank <- max_rank + 1
    }

    subtitle <- paste("n = ", n_model, ", p = ", p, ", k = ", k, ",
        beta_high = ", beta_high, sep="")

    plot <- ggplot(df_gg, aes(x=Model_Size, y=Weak_Signal_Selections,
        color=Label, shape=Label)) + geom_point(size=2.5, alpha=1) +
        xlab("Model Size") + ylab("Number of Weak Signal Features Selected") + 
        labs(subtitle=subtitle) + scale_x_continuous(breaks=seq(2, max_rank,
            by=2))

    if(!legend){
        plot <- plot + theme(legend.position="none")
    }

    return(plot)

}



createLossesPlot2 <- function(losses, names, j_choices_mat, n_model, p, k,
    beta_high, J=NA, legend=TRUE, xaxis="No. Fitted Coefficients",
    cluster_count_mat=NA){
    n_methods <- length(names)
    n_sims <- ncol(losses)

    # If J is not provided, figure out on own
    if(length(J) > 1){
        stop("length(J) > 1")
    }
    if(is.na(J)){
        J <- (nrow(losses) - 1)/n_methods
    }
    if(J != round(J)){
        print(nrow(losses))
        print(n_methods)
        print(J)
        stop("J != round(J)")
    }
    if(any(losses < 0)){
        stop("any(losses < 0)")
    }
    if(n_sims != ncol(j_choices_mat) | nrow(losses) != nrow(j_choices_mat) + 1){
        print(dim(losses))
        print(dim(j_choices_mat))
        stop("ncol(losses) != ncol(j_choices_mat) | nrow(losses) != nrow(j_choices_mat) + 1")
    }

    if(xaxis=="No. Selected Features"){
        if(n_sims != ncol(cluster_count_mat) | nrow(losses) != nrow(cluster_count_mat) + 1){
            print(dim(losses))
            print(dim(cluster_count_mat))
            stop("ncol(losses) != ncol(cluster_count_mat) | nrow(losses) != nrow(cluster_count_mat) + 1")
        }
        if(any(cluster_count_mat[j_choices_mat] == 0)){
            stop("any(cluster_count_mat[j_choices_mat] == 0)")
        }
    }

    if(xaxis=="No. Selected Features" & any(is.na(cluster_count_mat))){
        stop("If xaxis=No. Selected Features, must provide cluster_count_mat")
    }

    if(xaxis=="No. Fitted Coefficients"){

        if(n_sims > 1){
            row_inds <- which(rowSums(j_choices_mat) > 0)
            # if(!identical(which(rowSums(cluster_count_mat) > 0)), row_inds){
            #     stop("!identical(which(rowSums(cluster_count_mat) > 0)), row_inds")
            # }
            n_losses <- sum(rowSums(j_choices_mat) > 0)
            # if(!identical(sum(rowSums(cluster_count_mat) > 0), n_losses){
            #     stop("!identical(sum(rowSums(cluster_count_mat) > 0), n_losses")
            # }
            losses_vec <- numeric(n_losses)
            j_choices <- which(rowSums(j_choices_mat[1:J, ]) > 0)
            if(length(j_choices)*n_methods != n_losses){
                stop("length(j_choices)*n_methods + 1 != n_losses")
            }
            # if(!identical(which(rowSums(cluster_count_mat[1:J, ]) > 0),
            #     j_choices){
            #     stop("!identical(which(rowSums(cluster_count_mat[1:J, ]) > 0), j_choices")
            # }
            if(length(row_inds) != n_losses){
                stop("length(row_inds) != n_losses")
            }


            for(i in 1:n_losses){
                cols_i <- j_choices_mat[row_inds[i], ]
                losses_vec[i] <- mean(losses[row_inds[i], cols_i])
            }
        } else{
            row_inds <- which(j_choices_mat)
            n_losses <- sum(j_choices_mat)
            losses_vec <- numeric(n_losses)
            j_choices <- which(j_choices_mat[1:J, 1])
            if(length(j_choices)*n_methods != n_losses){
                stop("length(j_choices)*n_methods + 1 != n_losses")
            }
            if(length(row_inds) != n_losses){
                stop("length(row_inds) != n_losses")
            }

            losses_vec <- losses[row_inds, 1]
        }

        df_gg <- data.frame(rep(j_choices, n_methods), losses_vec,
            rep(names, each=length(j_choices)))

    } else if(xaxis=="No. Selected Features"){
        j_choices <- integer()
        losses_vec <- numeric()
        labels <- character()
        for(k in 1:n_methods){
            losses_vec_k <- as.vector(losses[((k - 1)*J + 1):(k*J), ])
            # print("summary(losses_vec_k:")
            # print(summary(losses_vec_k))
            # Get all model sizes for this method
            model_sizes <- as.vector(cluster_count_mat[((k - 1)*J + 1):(k*J), ])

            if(length(losses_vec_k) != length(model_sizes)){
                stop("length(losses_vec_k) != length(model_sizes)")
            }

            losses_vec_k <- losses_vec_k[model_sizes != 0]
            model_sizes <- model_sizes[model_sizes != 0]

            if(length(losses_vec_k) != length(model_sizes)){
                stop("length(losses_vec_k) != length(model_sizes)")
            }

            unique_model_sizes <- sort(unique(model_sizes))

            # loss_size_df <- as.data.frame(model_sizes, losses_vec_k)
            # colnames(loss_size_df) <- c("model_sizes", "losses_vec_k")
            # loss_size_df <- loss_size_df[!duplicated(loss_size_df$model_sizes, )]

            # sort(unique(

            
            n_sizes_k <- length(unique_model_sizes)
            # print("n_sizes_k:")
            # print(n_sizes_k)
            # print("unique_model_sizes:")
            # print(unique_model_sizes)



            j_choices <- c(j_choices, unique_model_sizes)
            labels <- c(labels, rep(names[k], n_sizes_k))
            losses_means_k <- numeric(n_sizes_k)
            for(j in 1:n_sizes_k){
                size_j <- unique_model_sizes[j]
                if(sum(model_sizes == size_j) == 0){
                    print("unique(model_sizes):")
                    print(unique(model_sizes))
                    print("size_j:")
                    print(size_j)
                    stop("sum(model_sizes == size_j) == 0")
                }
                # print("str(model_sizes == size_j):")
                # print(str(model_sizes == size_j))
                # print("losses_vec_k[model_sizes == size_j]:")
                # print(losses_vec_k[model_sizes == size_j])
                # print("str(losses_vec_k[model_sizes == size_j]):")
                # print(str(losses_vec_k[model_sizes == size_j]))
                # print("summary(losses_vec_k[model_sizes == size_j]):")
                # print(summary(losses_vec_k[model_sizes == size_j]))
                losses_means_k[j] <- mean(losses_vec_k[model_sizes == size_j],
                    na.rm=TRUE)

                if(is.na(losses_means_k[j])){
                    print("k:")
                    print(k)
                    print("losses_means_k[model_sizes == size_j]:")
                    print(losses_means_k[model_sizes == size_j])
                    stop("is.na(losses_means_k[j])")
                }
            }
            losses_vec <- c(losses_vec, losses_means_k)
            if(length(j_choices) != length(losses_vec) | length(losses_vec) != length(labels)){
                stop("length(j_choices) != length(losses_vec) | length(losses_vec) != length(labels)")
            }
        }
        df_gg <- data.frame(j_choices, losses_vec, labels)

        # print("str(df_gg):")
        # print(str(df_gg))
        # print("head(df_gg):")
        # print(head(df_gg))
    }

    

    true_model_loss <- mean(losses[nrow(losses), ])


    colnames(df_gg) <- c("Model_Size", "MSE", "Label")

    # Make labels display-friendly
    df_gg <- makeDfLabelsDisplay(df_gg)

    # Get max ranking value for integer labels on horizontal axis
    max_rank <- max(df_gg$Model_Size)
    if((max_rank %% 2) != 0){
        max_rank <- max_rank + 1
    }

    subtitle <- paste("n = ", n_model, ", p = ", p, ", k = ", k, ", beta_high = ",
        beta_high, sep="")

    plot <- ggplot(df_gg, aes(x=Model_Size, y=MSE, color=Label,
        shape=Label)) + geom_point(size=2.5, alpha=1) +
        geom_hline(yintercept=true_model_loss
            # , color="red"
            , linetype = "dashed"
            ) + xlab(xaxis) + labs(subtitle=subtitle) +
        scale_x_continuous(breaks=seq(2, max_rank, by=2))

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
    # to_avg: must be provided if average=TRUE. A logical vector of the same
    # length as selected. If feature j is in a cluster, the jth entry of 
    # to_avg will be TRUE. 
    #
    # avg_feats: must be provided if average=TRUE. A list of the same length as
    # selected. If feature j is in a
    # cluster, the jth entry of avg_feats will contain the features from the
    # cluster of which feature j is a member.
    #
    # weights: must be provided if average=TRUE. A list of the same length as
    # selected. If feature j is in a cluster, the jth entry of weights will
    # be the weights to use (in the same order as the jth entry of
    # avg_feats).

    p_selected <- length(selected)
    if(!average){
        X_df <- X[, selected]
    } else{
        X_df <- matrix(0, nrow(X), p_selected)
        colnames(X_df) <- character(p_selected)
        for(j in 1:p_selected){
            if(to_avg[j]){
                # Features to average
                feats_to_avg_j <- avg_feats[[j]]
                weights_j <- weights[[j]]
                if(abs(sum(weights_j) - 1) > 10^(-6)){
                    stop("abs(sum(weights_j) - 1) > 10^(-6)")
                }
                if(length(feats_to_avg_j) != length(weights_j)){
                    stop("length(feats_to_avg_j) != length(weights_j)")
                }
                # print("feats_to_avg_j:")
                # print(feats_to_avg_j)
                # print("weights_j:")
                # print(weights_j)
                # print("str(X[, feats_to_avg_j]):")
                # print(str(X[, feats_to_avg_j]))
                X_df[, j] <- X[, feats_to_avg_j] %*% weights_j
                if(p_selected > 1){
                    # print("str(X_df):")
                    # print(str(X_df))
                    # print(paste("X", selected[j], sep=""))
                    # print("colnames(X_df):")
                    # print(colnames(X_df))
                    # print("j:")
                    # print(j)
                    colnames(X_df)[j] <- paste("X", selected[j], sep="")
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

get_loss <- function(xtrain, ytrain, xtest, ytest, selected, average=FALSE,
    to_avg=NA, avg_feats=NA, weights=NA, digits=8) {
    # Inputs:
    #
    # selected: a vector of length j containing a selected set.
    #
    # average: if TRUE, will do an averaging of cluster members.
    #
    # to_avg: must be provided if average=TRUE. A logical vector of the same
    # length as selected. If feature j is in a cluster, the jth entry of 
    # to_avg will be TRUE. 
    #
    # avg_feats: must be provided if average=TRUE. A list of the same length as
    # selected. If feature j is in a
    # cluster, the jth entry of avg_feats will contain the features from the
    # cluster of which feature j is a member.
    #
    # weights: must be provided if average=TRUE. A list of the same length as
    # selected. If feature j is in a cluster, the jth entry of weights will
    # be the weights to use (in the same order as the jth entry of
    # avg_feats).
    #
    p_selected <- length(selected)
    if(p_selected == 0){
        stop("length(selected) == 0")
    }
    n_train <- nrow(xtrain)
    n_test <- nrow(xtest)
    p <- ncol(xtrain)
    if(length(ytrain) != n_train){
        stop("length(ytrain) != n_train")
    }
    if(length(ytest) != n_test){
        stop("length(ytest) != n_test")
    }
    if(p != ncol(xtest)){
        stop("p != ncol(xtest)")
    }
    if(max(selected) > p){
        stop("max(selected) > p")
    }
    if(p_selected > p){
        stop("length(selected) > p")
    }
    if(average){
        if(all(is.na(to_avg))){
            stop("If average==TRUE, must provide to_avg to get_loss")
        }
        if(any(to_avg)){
            if(all(is.na(avg_feats)) | all(is.na(weights))){
                stop("If average==TRUE, must provide avg_feats and weights to get_loss")
            }
        }
        if(length(to_avg) > p){
            stop("length(to_avg) > ncol(xtrain)")
        }
        if(!all.equal(lengths(avg_feats), lengths(weights))){
            stop("!all.equal(lengths(avg_feats), lengths(weights))")
        }
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
    loss <- Metrics::mse(actual=y_test, predicted=predict(model,
        newdata=data_test))
    return(round(loss, digits))
}

num_false <- function(selected){
    if(length(selected) == 0){
        setwd(wd)
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
        setwd(wd)
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
    if(method_index == 0){
        stop("method_index == 0")
    }
    return(method_index)
}


extractSelectedSets <- function(out, method, j_choices){
    # Extracts selected sets given an @out object from simulator and
    # the name of the method

    # List of selected sets to return
    selected_sets <- list()
    if(length(j_choices) == 0){
        warning("j_choices is empty")
        return(list(selected_sets, j_choices))
    }

    # J's to elminate from j_choices
    j_choices_to_eliminate <- integer()

    # Results from this iteration

    
    if(method %in% c("lasso", "lasso_random", "lasso_proto")){
        lasso_sets <- unique(out$selected)
        for(j in j_choices){
            # Lasso selected set of size j
            lasso_sets_j <- lasso_sets[lapply(lasso_sets, length) == j]
            if(length(lasso_sets_j) == 0){
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
            selected_sets[[j]] <- lasso_sets_j[[1]]
        }
    } else if(method == "BRVZ_avg_unwt"){
        for(j in j_choices){
            j_choices <- j_choices[j_choices <= max(vapply(out$selected_sets,
                length, integer(1)))]
            # Also need: to_avg, avg_feats, weights
            to_avg_list <- list()
            avg_feats_list <- list()
            weights_list <- list()
            for(j in j_choices){
                if(length(out$selected_sets[[j]]) > 0){
                    # The first j selected features
                    selected_sets[[j]] <- out$selected_sets[[j]]
                    # to_avg_list: A list whose elements are logical vectors of length j.
                    # If feature k in 1:j is in a 
                    # cluster, the kth entry of the jth vector in to_avg_list will be
                    # TRUE. 
                    to_avg_list[[j]] <- out$to_avg_list[[j]]
                    # avg_feats_list: A list whose elements are lists of length j. If
                    # feature k in 
                    # 1:j is in a cluster, the kth entry of the jth vector of avg_feats_list
                    # will contain the features from the
                    # cluster of which feature k is a member.
                    avg_feats_list[[j]] <- out$avg_feats_list[[j]]
                    # weights_list: A list whose elements are lists of length j. If
                    # feature k in 
                    # 1:j is in a cluster, the kth entry of the jth vector of
                    # weights_list will
                    # be the weights to use (in the same order as the jth entry of
                    # avg_feats_list).
                    weights_list[[j]] <- out$weights_list[[j]]
                } else{
                    j_choices_to_eliminate <- c(j_choices_to_eliminate, j)
                }
            }
        }
    } else if(method %in% c("lassoSS_phat", "lassoMB_phat",
        "SS_SS_random")){
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

            if(stability_selected_top_j_tau >= decreasing_taus[j] |
                stability_selected_top_j_tau <= decreasing_taus[j + 1]){
                stop("stability_selected_top_j_tau >= decreasing_taus[j] | stability_selected_top_j_tau <= decreasing_taus[j + 1]")
            }

            selected_sets[[j]] <- which(taus >
                stability_selected_top_j_tau)

        }
    } else if(method %in% c("subspace_sel_grass", "subspace_ss_max_cancor",
        "subspace_ss_min_cancor")){
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
        # feaatures in decreasing order of weighted selection proportion
        # by default
        # Selected sets
        selected <- out$selected
        # print("selected:")
        # print(selected)
        j_choices <- j_choices[j_choices <= length(selected)]
        for(j in j_choices){
            selected_sets[[j]] <- selected[1:j]
        }
        # print("selected_sets:")
        # print(selected_sets)
    } else if(method %in% c("SS_GSS_random_avg", "SS_GSS_random_avg_unwt")){
        # Selected sets
        selected <- out$selected
        j_choices <- j_choices[j_choices <= length(selected)]
        # Also need: to_avg, avg_feats, weights
        to_avg_list <- list()
        avg_feats_list <- list()
        weights_list <- list()
        for(j in j_choices){
            # The first j selected features
            selected_sets[[j]] <- selected[1:j]
            # to_avg_list: A list whose elements are logical vectors of length j.
            # If feature k in 1:j is in a 
            # cluster, the kth entry of the jth vector in to_avg_list will be
            # TRUE. 
            to_avg_list[[j]] <- out$to_avg[1:j]
            # avg_feats_list: A list whose elements are lists of length j. If
            # feature k in 
            # 1:j is in a cluster, the kth entry of the jth vector of avg_feats_list
            # will contain the features from the
            # cluster of which feature k is a member.
            avg_feats_list[[j]] <- out$avg_feats[1:j]
            # weights_list: A list whose elements are lists of length j. If
            # feature k in 
            # 1:j is in a cluster, the kth entry of the jth vector of
            # weights_list will
            # be the weights to use (in the same order as the jth entry of
            # avg_feats_list).
            weights_list[[j]] <- out$weights[1:j]
        }
        
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

    if(!(method %in% c("SS_GSS_random_avg", "SS_GSS_random_avg_unwt",
        "BRVZ_avg_unwt"))){
        return(list(selected_sets, j_choices))
    }

    # Outputs:
    #
    # to_avg_list: A list whose elements are logical vectors of length j.
    # If feature k in 1:j is in a 
    # cluster, the kth entry of the jth vector in to_avg_list will be
    # TRUE. 
    #
    # avg_feats_list: A list whose elements are lists of length j. If
    # feature k in 
    # 1:j is in a cluster, the kth entry of the jth vector of avg_feats_list
    # will contain the features from the
    # cluster of which feature k is a member.
    #
    # weights_list: A list whose elements are lists of length j. If
    # feature k in 
    # 1:j is in a cluster, the kth entry of the jth vector of
    # weights_list will
    # be the weights to use (in the same order as the jth entry of
    # avg_feats_list).

    
    

    return(list(selected_sets, j_choices, to_avg_list, avg_feats_list,
        weights_list))
}


calcHamStab <- function(sel_mat, n_sims){
    require(e1071)
    if(nrow(sel_mat) != n_sims){
        stop("nrow(sel_mat) != n_sims")
    }
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
    if(ncol(sel_mat_nonzero) != p){
        stop("ncol(sel_mat_nonzero) != p")
    }
    m <- nrow(sel_mat_nonzero)
    if(m == 1){
        stop("m == 1")
    }
    if(m <= n_sims/2){
        print("Warning here")
        warning("m <= n_sims/2")
    }
    ham_dist_mat <- e1071::hamming.distance(sel_mat_nonzero)
    return(sum(ham_dist_mat[lower.tri(ham_dist_mat)])/(m*(m-1)/2))
}

calcNogueiraStab <- function(sel_mat, n_sims, calc_errors=FALSE, conf=0.95,
    coarseness=1, cluster_count=FALSE, n_clust_feats=0){
    if(nrow(sel_mat) != n_sims*coarseness){
        stop("nrow(sel_mat) != n_sims*coarseness")
    }

    if(cluster_count & (n_clust_feats==0)){
        stop("Must provide n_clust_feats if cluster_count==TRUE")
    }


    p <- ncol(sel_mat)
    if(all(sel_mat == 0)){
        return(NA)
    }
    # Eliminate rows of all zeroes (no selections)
    all_zero_rows <- apply(sel_mat, 1, function(x){all(x == 0)})
    M <- sum(!all_zero_rows)
    if(M == 1){
        return(NA)
    }

    if(M <= n_sims*coarseness/10){
        print("Warning here")
        warning("sum(!all_zero_rows) <= n_sims*coarseness/10")
        return(NA)
    }
    sel_mat_nonzero <- sel_mat[!all_zero_rows, ]
    if(ncol(sel_mat_nonzero) != p){
        stop("ncol(sel_mat_nonzero) != p")
    }

    if(cluster_count){
        # If modifying criterion for clusters: for each row, count how
        # many clustered features were selected
        n_clust_selecs <- rowSums(sel_mat_nonzero[, 1:n_clust_feats])
        sel_mat_nonzero[, 1:n_clust_feats] <- 0
        for(i in 1:nrow(sel_mat_nonzero)){
            # For the ith row, fill the first n_clust_selecs[i] cluster
            # features with ones (i.e., don't distinguish amongst clustered
            # features)
            sel_mat_nonzero[i, 1:n_clust_selecs[i]] <- 1
        }
    }

    # Calculate s_f_squared values
    s_f_squared <- apply(sel_mat_nonzero, 2, var)
    k <- rowSums(sel_mat_nonzero)
    k_bar <- mean(k)

    stat <- 1 - mean(s_f_squared)/(k_bar/p*(1 - k_bar/p))

    if(calc_errors){
        alpha <- 1 - conf
        p_hat <- colMeans(sel_mat_nonzero,)
        # Calculate phi_hat_i's
        phi_hat_i <- numeric(M)
        for(i in 1:M){
            phi_hat_i[i] <- 1/(k_bar/p*(1 - k_bar/p))*
                (mean(sel_mat_nonzero[i, ]*p_hat) - k[i]*k_bar/p^2 + stat/2*(
                    2*k_bar*k[i]/p^2 - k[i]/p - k_bar/p + 1))
        }
        if(any(phi_hat_i == 0)){
            stop("any(phi_hat_i == 0)")
        }
        v <- 4/M^2*sum((phi_hat_i - mean(phi_hat_i))^2)
        if(v == 0){
            return(rep(NA, 3))
        }
        if(v < 0){
            stop("v < 0")
        }
        # Calculate errors for error bars
        margin <- qnorm(1 - alpha/2)*sqrt(v)
        return(c(stat, stat - margin, stat + margin))
    }

    return(stat)

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

    colnames(df_gg) <- c("Model_Size", "HammingInstability", "Label")

    # Make labels display-friendly
    df_gg <- makeDfLabelsDisplay(df_gg)

    # Get max ranking value for integer labels on horizontal axis
    max_rank <- max(df_gg$Model_Size)
    if((max_rank %% 2) != 0){
        max_rank <- max_rank + 1
    }

    subtitle <- paste("n = ", n_model, ", p = ", p, ", k = ", k, ",
        beta_high = ", beta_high, sep="")

    plot <- ggplot(df_gg, aes(x=Model_Size, y=HammingInstability, color=Label,
        shape=Label)) + geom_point(size=2.5, alpha=1) +
        xlab("Model Size") + ylab("Hamming Instability") +
        labs(subtitle=subtitle) +
        scale_x_continuous(breaks=seq(2, max_rank, by=2))

    if(!legend){
        plot <- plot + theme(legend.position="none")
    }

    return(plot)
}

createNogueiraStabPlot <- function(stab_mets, n_model, p, k,
    beta_high, legend=TRUE, plot_errors=FALSE, lowers=NA, uppers=NA,
    cluster_count=FALSE){
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
    if(plot_errors){
        lowers <- lowers[!all_na_rows, ]
        uppers <- uppers[!all_na_rows, ]
    }

    df_gg <- data.frame(rep(j_choices, n_methods), as.vector(stab_mets),
        rep(names, each=length(j_choices)))

    colnames(df_gg) <- c("Model_Size", "NogueiraStability", "Label")

    # Make labels display-friendly
    df_gg <- makeDfLabelsDisplay(df_gg)

    # Get max ranking value for integer labels on horizontal axis
    max_rank <- max(df_gg$Model_Size)
    if((max_rank %% 2) != 0){
        max_rank <- max_rank + 1
    }

    subtitle <- paste("n = ", n_model, ", p = ", p, ", k = ", k, ",
        beta_high = ", beta_high, sep="")

    if(cluster_count){
        subtitle <- paste("Cluster counting,", subtitle)
    }

    if(plot_errors){
        # Add in error bounds
        lower_bounds <- as.vector(lowers)
        upper_bounds <- as.vector(uppers)
        if(length(lower_bounds) != nrow(df_gg) |
            length(upper_bounds) != nrow(df_gg)){
            stop("length(lower_bounds) != nrow(df_gg) | length(upper_bounds) != nrow(df_gg)")
        }
        df_gg <- data.frame(df_gg, lower_bounds, upper_bounds)
        colnames(df_gg)[(ncol(df_gg) - 1):ncol(df_gg)] <- c("Lower", "Upper")
    }

    plot <- ggplot(df_gg, aes(x=Model_Size, y=NogueiraStability, color=Label,
        shape=Label)) + geom_point(size=2.5, alpha=1) +
        xlab("Model Size") + ylab("Nogueira Stability") + labs(subtitle=subtitle) +
        scale_x_continuous(breaks=seq(2, max_rank, by=2))

    if(plot_errors){
        plot <- plot + geom_errorbar(aes(ymin = Lower, ymax = Upper),
            width = 0.5)
    }

    if(!legend){
        plot <- plot + theme(legend.position="none")
    }



    return(plot)
}







