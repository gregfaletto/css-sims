# setwd("/Users/gregfaletto/Google Drive/Data Science/LaTeX/Paper skeleton 2020/Real data examples")



library(ggplot2)
library(Metrics)
# library(stabs)

# Directory where stabsel_copy.R is storeed
dir_stabsel <- "/Users/gregfaletto/Google Drive/Data Science/R/USC/Stability Selection/toy_example/New method"

wd <- getwd()
setwd(dir_stabsel)
source(file="stabsel_copy.R")
setwd(wd)

genResponse <- function(data, snr, beta_high, latent_feat, cluster_feats,
    unclustered_feats=NA){
    p_clustered <- length(cluster_feats)
    # p_unclustered <- length(unclustered_feats)
    n <- nrow(data)
    if(length(latent_feat) != n){
        stop("length(latent_feat) != nrow(data)")
    }
    if(!all(cluster_feats %in% colnames(data))){
        stop("!all(cluster_feats %in% colnames(data))")
    }
    # if(!all(unclustered_feats %in% colnames(data))){
    #     stop("!all(unclustered_feats %in% colnames(data))")
    # }
    # Generate response
    
    
    mu <- beta_high*latent_feat +
        data[, cluster_feats] %*% rep(1, p_clustered) 
        # + data[, unclustered_feats] %*% rep(1, p_unclustered)
    sd <- sqrt(sum(mu^2) / (n * snr))
    y <- mu + sd*rnorm(n)
    return(list(mu, y))
}

genIndices <- function(n, n_selec, n_train){
    n_test <- n - n_train - n_selec
    if(n_test <= 0){
        stop("n_test <= 0")
    }
    # Randomly divide into selection, training, and test sets
    selec_indices <- sample(1:n, n_selec)
    train_indices <- sample(setdiff(1:n, selec_indices), n_train)
    if(length(intersect(selec_indices, train_indices)) > 0){
        stop("length(intersect(selec_indices, train_indices)) > 0")
    }
    test_indices <- setdiff(1:n, c(selec_indices, train_indices))
    if(length(test_indices) != n_test){
        stop("length(test_indices) != n_test")
    }
    return(list(selec_indices, train_indices, test_indices))
}

genIndicesNoTest <- function(n, n_selec){
    n_train <- n - n_selec
    if(n_train <= 0){
        stop("n_train <= 0")
    }
    # Randomly divide into selection, training, and test sets
    selec_indices <- sample(1:n, n_selec)
    train_indices <- setdiff(1:n, selec_indices)
    if(length(train_indices) != n_train){
        stop("length(train_indices) != n_train")
    }
    return(list(selec_indices, train_indices))
}

# createLossesPlot3 <- function(losses, j_choices, coarseness=1){
#     if(any(losses[complete.cases(losses), ] < 0)){
#         stop("any(losses[complete.cases(losses), ] < 0)")
#     }
#     if(!all(j_choices %in% 1:nrow(losses))){
#         stop("!all(j_choices %in% 1:nrow(losses))")
#     }

#     losses <- losses[j_choices, ]
#     if(!all(complete.cases(losses))){
#         stop("!all(complete.cases(losses))")
#     }

#     losses_vec <- c(losses[, 1], losses[, 2], losses[, 3])

#     df_gg <- data.frame(rep(j_choices, 3), losses_vec, rep(colnames(losses),
#     	each=length(j_choices)))

#     colnames(df_gg) <- c("Model_Size", "MSE", "Method")

#     plot <- ggplot(df_gg, aes(x=Model_Size, y=MSE, color=Method,
#         shape=Method)) + geom_point(size=2.5, alpha=1)+ xlab("Model Size")

#     return(plot)
# }

makeDfMethodsDisplay <- function(df){
    # Takes in data.frame with computer-friendly labels and returns 
    # data.frame with display-friendly labels
    if(!("Method" %in% colnames(df))){
        stop("!(Method %in% colnames(df))")
    }
    if(ncol(df) != 3){
        stop("ncol(df) != 3 (before modifications)")
    }
    labels <- as.factor(nameMap(as.character(df$Method)))
    labels <- droplevels(labels)
    df$Method <- labels
    if(!all(levels(df$Method) %in% c("Lasso", "Stability Selection",
        "Sparse CSS"))){
        stop("!all(levels(df$Method) %in% c(Lasso, Stability Selection,
            Sparse CSS))")
    }
    if(ncol(df) != 3){
        stop("ncol(df) != 3 (after modifications)")
    }
    return(df)
}

makeDfMethodsDisplay3 <- function(df){
    # Takes in data.frame with computer-friendly labels and returns 
    # data.frame with display-friendly labels
    if(!("Method" %in% colnames(df))){
        stop("!(Method %in% colnames(df))")
    }
    if(ncol(df) != 3){
        stop("ncol(df) != 3 (before modifications)")
    }
    labels <- as.factor(nameMap(as.character(df$Method)))
    labels <- droplevels(labels)
    df$Method <- labels
    if(!all(levels(df$Method) %in% c("Lasso", "Stability Selection",
        "Sparse CSS", "Weighted Averaged CSS",
        "Simple Averaged CSS"))){
        stop("!all(levels(df$Method) %in% c(Lasso, Stability Selection,
            Sparse CSS), Weighted Averaged CSS, Simple Averaged CSS)")
    }
    if(ncol(df) != 3){
        stop("ncol(df) != 3 (after modifications)")
    }
    return(df)
}


createLossesPlot4 <- function(losses, names, j_choices_mat, J=NA){
    n_methods <- length(names)
    n_sims <- ncol(losses)

    # If J is not provided, figure out on own
    if(length(J) > 1){
        stop("length(J) > 1")
    }
    if(is.na(J)){
        J <- (nrow(losses))/n_methods
    }
    if(J != round(J)){
        print(nrow(losses))
        print(n_methods)
        print(J)
        stop("J != round(J)")
    }

    if(n_sims != ncol(j_choices_mat) | nrow(losses) != nrow(j_choices_mat)){
        print(dim(losses))
        print(dim(j_choices_mat))
        stop("ncol(losses) != ncol(j_choices_mat) | nrow(losses) != nrow(j_choices_mat)")
    }

    if(n_sims > 1){
        row_inds <- which(rowSums(j_choices_mat) > 0)
        n_losses <- sum(rowSums(j_choices_mat) > 0)
        losses_vec <- numeric(n_losses)
        j_choices <- which(rowSums(j_choices_mat[1:J, ]) > 0)
        if(length(j_choices)*n_methods != n_losses){
            stop("length(j_choices)*n_methods + 1 != n_losses")
        }
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
            print(J)
            print(j_choices)
            print(j_choices_mat)
            print(dim(j_choices_mat))
            print(length(j_choices))
            print(n_methods)
            print(n_losses)
            stop("length(j_choices)*n_methods != n_losses")
        }
        if(length(row_inds) != n_losses){
            stop("length(row_inds) != n_losses")
        }

        losses_vec <- losses[row_inds, 1]
    }

    df_gg <- data.frame(rep(j_choices, n_methods), losses_vec,
        rep(names, each=length(j_choices)))

    colnames(df_gg) <- c("Model_Size", "MSE", "Method")

    plot <- ggplot(df_gg, aes(x=Model_Size, y=MSE, color=Method,
        shape=Method)) + geom_point(size=2.5, alpha=1) +
        xlab("Model Size")

    return(plot)
}

createProxiesPlot4 <- function(proxy_counts, names, j_choices_mat, J=NA){
    n_methods <- length(names)
    n_sims <- ncol(proxy_counts)

    # If J is not provided, figure out on own
    if(length(J) > 1){
        stop("length(J) > 1")
    }
    if(is.na(J)){
        J <- (nrow(proxy_counts))/n_methods
    }
    if(J != round(J)){
        print(nrow(proxy_counts))
        print(n_methods)
        print(J)
        stop("J != round(J)")
    }

    if(n_sims != ncol(j_choices_mat) | nrow(proxy_counts) != nrow(j_choices_mat)){
        print(dim(proxy_counts))
        print(dim(j_choices_mat))
        stop("ncol(proxy_counts) != ncol(j_choices_mat) | nrow(proxy_counts) != nrow(j_choices_mat)")
    }

    if(n_sims > 1){
        row_inds <- which(rowSums(j_choices_mat) > 0)
        n_proxy_counts <- sum(rowSums(j_choices_mat) > 0)
        proxy_counts_vec <- numeric(n_proxy_counts)
        j_choices <- which(rowSums(j_choices_mat[1:J, ]) > 0)
        if(length(j_choices)*n_methods != n_proxy_counts){
            stop("length(j_choices)*n_methods + 1 != n_proxy_counts")
        }
        if(length(row_inds) != n_proxy_counts){
            stop("length(row_inds) != n_proxy_counts")
        }


        for(i in 1:n_proxy_counts){
            cols_i <- j_choices_mat[row_inds[i], ]
            proxy_counts_vec[i] <- mean(proxy_counts[row_inds[i], cols_i])
        }
    } else{
        row_inds <- which(j_choices_mat)
        n_proxy_counts <- sum(j_choices_mat)
        proxy_counts_vec <- numeric(n_proxy_counts)
        j_choices <- which(j_choices_mat[1:J, 1])
        if(length(j_choices)*n_methods != n_proxy_counts){
            print(J)
            print(j_choices)
            print(j_choices_mat)
            print(dim(j_choices_mat))
            print(length(j_choices))
            print(n_methods)
            print(n_proxy_counts)
            stop("length(j_choices)*n_methods != n_proxy_counts")
        }
        if(length(row_inds) != n_proxy_counts){
            stop("length(row_inds) != n_proxy_counts")
        }

        proxy_counts_vec <- proxy_counts[row_inds, 1]
    }

    df_gg <- data.frame(rep(j_choices, n_methods), proxy_counts_vec,
        rep(names, each=length(j_choices)))

    colnames(df_gg) <- c("Model_Size", "Strong_Proxies", "Method")

    plot <- ggplot(df_gg, aes(x=Model_Size, y=Strong_Proxies, color=Method,
        shape=Method)) + geom_point(size=2.5, alpha=1) +
        xlab("Model Size")

    return(plot)
}

createAllProxiesPlot4 <- function(all_proxy_counts, names, j_choices_mat, J=NA){
    n_methods <- length(names)
    n_sims <- ncol(all_proxy_counts)

    # If J is not provided, figure out on own
    if(length(J) > 1){
        stop("length(J) > 1")
    }
    if(is.na(J)){
        J <- (nrow(all_proxy_counts))/n_methods
    }
    if(J != round(J)){
        print(nrow(all_proxy_counts))
        print(n_methods)
        print(J)
        stop("J != round(J)")
    }

    if(n_sims != ncol(j_choices_mat) | nrow(all_proxy_counts) != nrow(j_choices_mat)){
        print(dim(all_proxy_counts))
        print(dim(j_choices_mat))
        stop("ncol(all_proxy_counts) != ncol(j_choices_mat) | nrow(all_proxy_counts) != nrow(j_choices_mat)")
    }

    if(n_sims > 1){
        row_inds <- which(rowSums(j_choices_mat) > 0)
        n_all_proxy_counts <- sum(rowSums(j_choices_mat) > 0)
        all_proxy_counts_vec <- numeric(n_all_proxy_counts)
        j_choices <- which(rowSums(j_choices_mat[1:J, ]) > 0)
        if(length(j_choices)*n_methods != n_all_proxy_counts){
            stop("length(j_choices)*n_methods + 1 != n_all_proxy_counts")
        }
        if(length(row_inds) != n_all_proxy_counts){
            stop("length(row_inds) != n_all_proxy_counts")
        }


        for(i in 1:n_all_proxy_counts){
            cols_i <- j_choices_mat[row_inds[i], ]
            all_proxy_counts_vec[i] <- mean(all_proxy_counts[row_inds[i], cols_i])
        }
    } else{
        row_inds <- which(j_choices_mat)
        n_all_proxy_counts <- sum(j_choices_mat)
        all_proxy_counts_vec <- numeric(n_all_proxy_counts)
        j_choices <- which(j_choices_mat[1:J, 1])
        if(length(j_choices)*n_methods != n_all_proxy_counts){
            print(J)
            print(j_choices)
            print(j_choices_mat)
            print(dim(j_choices_mat))
            print(length(j_choices))
            print(n_methods)
            print(n_all_proxy_counts)
            stop("length(j_choices)*n_methods != n_all_proxy_counts")
        }
        if(length(row_inds) != n_all_proxy_counts){
            stop("length(row_inds) != n_all_proxy_counts")
        }

        all_proxy_counts_vec <- all_proxy_counts[row_inds, 1]
    }

    df_gg <- data.frame(rep(j_choices, n_methods), all_proxy_counts_vec,
        rep(names, each=length(j_choices)))

    colnames(df_gg) <- c("Model_Size", "All_Proxies", "Method")

    plot <- ggplot(df_gg, aes(x=Model_Size, y=All_Proxies, color=Method,
        shape=Method)) + geom_point(size=2.5, alpha=1) +
        xlab("Model Size")

    return(plot)
}

get_loss <- function(xtrain, ytrain, xtest, mutest, selected, digits=8){
    if(length(selected) == 0){
        stop("length(selected) == 0")
    }

    data_train <- data.frame(ytrain, xtrain[, selected])
    colnames(data_train)[1] <- "y"

    # If only one feature was selected, then
    # x_train[, selected] is a vector, and the name of this column
    # of the data.frame will be weird. need to re-name the second column of
    # the data.frame in this case.
    if(length(selected) == 1){
        colnames(data_train)[2] <- paste("X", selected, sep="")
    }

    data_test <- data.frame(mutest, xtest[, selected])
    colnames(data_test)[1] <- "y"

    if(length(selected) == 1){
        colnames(data_test)[2] <- paste("X", selected, sep="")
    }
    if(any(colnames(data_test)[2:ncol(data_test)] !=
        colnames(data_train)[2:ncol(data_train)])){
        colnames(data_test)[2:ncol(data_test)] <-
            colnames(data_train)[2:ncol(data_train)]
    }

    model <- lm(y~., data_train)
    loss <- Metrics::mse(actual=mutest, predicted=predict(model,
        newdata=data_test))
    return(round(loss, digits))
}

get_loss_no_test <- function(xtrain, ytrain, mutrain, selected, digits=8){
    if(length(selected) == 0){
        stop("length(selected) == 0")
    }

    data_train <- data.frame(ytrain, xtrain[, selected])
    colnames(data_train)[1] <- "y"

    # If only one feature was selected, then
    # x_train[, selected] is a vector, and the name of this column
    # of the data.frame will be weird. need to re-name the second column of
    # the data.frame in this case.
    if(length(selected) == 1){
        colnames(data_train)[2] <- paste("X", selected, sep="")
    }

    model <- lm(y~., data_train)
    loss <- Metrics::mse(actual=mutrain, predicted=predict(model))
    return(round(loss, digits))
}

getLassoSelected <- function(X_train, y_train, j_choices){
	selected_sets <- list(max(j_choices))
	# Fit lasso
	lasso <- glmnet(X_train, y_train, nlambda=5000)
	lasso_sets <- unique(predict(lasso, type="nonzero"))
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
            j_choices <- setdiff(j_choices, j)
            next
        }
        # print(paste("Lasso selected set of size", j, ":"))
        # print(colnames(X_train)[lasso_sets_j[[1]]])
        selected_sets[[j]] <- lasso_sets_j[[1]]
    }
    return(list(selected_sets, j_choices))
}

getStabSelected <- function(X_train, y_train, j_choices, q, B){
	selected_sets <- list(max(j_choices))
	# Get stability selected sets
	fit <- stabs::stabsel(x=X_train, y=y_train, fitfun=glmnet.lasso,
		cutoff=0.7, q=q, B=B, sampling_type="SS")
	p_hat <- fit$feat_sel_props
    decreasing_taus <- sort(p_hat, decreasing=TRUE)
    # print("decreasing_taus (stability selection):")
    # print(decreasing_taus)
    j_max <- sum(decreasing_taus > 0.5)
    j_choices <- intersect(j_choices, 1:j_max)
    for(j in j_choices){

        # First, check if we need to skip this j. True if there is a tie in 
        # probabilities resulting in no stability selected set of size j
        # being defined.

        if(decreasing_taus[j] == decreasing_taus[j + 1]){
            # If there is a tie in probabilities, then there is no stability
            # selected set of size j defined. In that case, make a note of
            # this entry (will need to delete it later from both vectors) and
            #  move on to the next j.
            j_choices <- setdiff(j_choices, j)
            next
        }
        
        # Stability selected set of size j    

        stability_selected_top_j_tau <- (decreasing_taus[j] +
            decreasing_taus[j + 1])/2

        stability_selected_sets_j  <- which(p_hat >
            stability_selected_top_j_tau)

        selected_sets[[j]] <- stability_selected_sets_j
    }
    return(list(selected_sets, j_choices))
}

getGenStabSelected <- function(X_train, y_train, j_choices, R, q, B){
	selected_sets <- list(max(j_choices))
	# Get stability selected sets
	
	fit <- stabselGreg(x=X_train, y=y_train, fitfun=glmnet.lasso,
		cutoff=0.501, q=q, B=B, sampling_type="SS", papply = mclapply, R=R)

	selected <- fit$selected

	# print("selected (gen stab sel):")
	# print(colnames(X_train)[selected])

    j_choices <- j_choices[j_choices <= length(selected)]
    for(j in j_choices){
        selected_sets[[j]] <- selected[1:j]
    }

    return(list(selected_sets, j_choices))
}

dataResults <- function(X_selec, y_selec, X_train, y_train, X_test, mu_test, R,
	p_max){
	# Possible sizes of selected sets
    j_choices <- 1:p_max
	# Get selected sets
	lasso_results <- getLassoSelected(X_selec, y_selec, j_choices)
	lasso_selected <- lasso_results[[1]]
	j_choices <- lasso_results[[2]]
	stab_results <- getStabSelected(X_selec, y_selec, j_choices, q=p_max,
		B=1000)
	stab_selected <- stab_results[[1]]
	j_choices <- stab_results[[2]]
	gen_stab_results <- getGenStabSelected(X_selec, y_selec, j_choices, R=R,
		q=p_max, B=1000)
	gen_stab_selected <- gen_stab_results[[1]]
	j_choices <- gen_stab_results[[2]]

	# Fit OLS fits for lasso selected set, stability selection selected set,
    # and generalized stability selected set
    losses_mat <- matrix(NA, nrow=p_max, ncol=3)
    if(!all(j_choices %in% 1:p_max)){
        stop("!all(j_choices %in% 1:p_max)")
    }
    for(j in j_choices){
        # losses[j, k] <- get_loss(x_train, y_train, x_test,
        #     y_test, selected_sets[[k]][[j]])
        losses_mat[j, 1] <- get_loss(X_train, y_train,
            X_test, mu_test, lasso_selected[[j]])
        losses_mat[j, 2] <- get_loss(X_train, y_train,
            X_test, mu_test, stab_selected[[j]])
        losses_mat[j, 3] <- get_loss(X_train, y_train,
            X_test, mu_test, gen_stab_selected[[j]])
    }

    # losses_mat <- losses_mat[complete.cases(losses_mat), ]
    colnames(losses_mat) <- c("Lasso", "Stability Selection",
    	"Sparse CSS")

    if(!all(j_choices %in% 1:nrow(losses_mat))){
        stop("!all(j_choices %in% 1:nrow(losses_mat))")
    }

    plot_overall2 <- createLossesPlot3(losses_mat, j_choices)

    print(plot_overall2)
}

dataResultsLoop <- function(X_selec, y_selec, X_train, y_train, X_test, mu_test,
    R, p_max, losses=NA, j_choices=NA){
    # Current number of columns in losses (needs to be 1 more than this later)

    if(all(is.na(losses))){
        n_cols_losses <- 0
    } else{
        n_cols_losses <- ncol(losses)
    }

    if(n_cols_losses != 0){
        if (nrow(losses) != 3*p_max){
            print("p_max:")
            print(p_max)
            print("dim(losses):")
            print(dim(losses))
            stop("nrow(losses) != 3*p_max")
        }
    }

    if(all(is.na(j_choices))){
        n_cols_j_choices <- 0
    } else{
        n_cols_j_choices <- ncol(j_choices)
    }

    # Possible sizes of selected sets
    j_choices_range <- 1:p_max
    # Get selected sets
    lasso_results <- getLassoSelected(X_selec, y_selec, j_choices_range)
    lasso_selected <- lasso_results[[1]]
    j_choices_range <- lasso_results[[2]]
    stab_results <- getStabSelected(X_selec, y_selec, j_choices_range, q=p_max,
        B=100)
    stab_selected <- stab_results[[1]]
    j_choices_range <- stab_results[[2]]
    gen_stab_results <- getGenStabSelected(X_selec, y_selec, j_choices_range, R=R,
        q=p_max, B=100)
    gen_stab_selected <- gen_stab_results[[1]]
    j_choices_range <- gen_stab_results[[2]]

    # Fit OLS fits for lasso selected set, stability selection selected set,
    # and generalized stability selected set
    losses_mat <- matrix(NA, nrow=3*p_max, ncol=1)
    if(!all(j_choices_range %in% 1:p_max)){
        stop("!all(j_choices_range %in% 1:p_max)")
    }
    for(j in j_choices_range){
        # losses[j, k] <- get_loss(x_train, y_train, x_test,
        #     y_test, selected_sets[[k]][[j]])
        losses_mat[j, 1] <- get_loss(X_train, y_train,
            X_test, mu_test, lasso_selected[[j]])
        losses_mat[j + p_max, 1] <- get_loss(X_train, y_train,
            X_test, mu_test, stab_selected[[j]])
        losses_mat[j + 2*p_max, 1] <- get_loss(X_train, y_train,
            X_test, mu_test, gen_stab_selected[[j]])
    }


    if(all(is.na(losses))){
        losses_ret <- losses_mat
    } else{
        losses_ret <- cbind(losses, losses_mat)
    }

    if(ncol(losses_ret) != n_cols_losses + 1){
        print("ncol(losses_ret):")
        print(ncol(losses_ret))
        print("n_cols_losses:")
        print(n_cols_losses)
        stop("ncol(losses_ret) != n_cols_losses + 1")
    }

    j_choices_mat <- matrix(FALSE, nrow=3*p_max, ncol=1)
    for(i in 1:3){
        j_choices_mat[j_choices_range + p_max*(i-1), 1] <- TRUE
    }
    j_choices_mat <- as.matrix(j_choices_mat)
    if(!all(j_choices_mat[1:p_max, ] ==
        j_choices_mat[(p_max + 1):(2*p_max), ]) |
        !all(j_choices_mat[1:p_max, ] ==
        j_choices_mat[(2*p_max + 1):(3*p_max), ])){
        print(j_choices_mat)
        stop("!all.equal(j_choices_mat[1:p_max, ],
        j_choices_mat[(p_max + 1):(2*p_max), ]) |
        !all.equal(j_choices_mat[1:p_max, ],
        j_choices_mat[(2*p_max + 1):(3*p_max), ])")
    }

    if(all(is.na(losses))){
        losses_ret <- losses_mat
    } else{
        losses_ret <- cbind(losses, losses_mat)
    }


    if(n_cols_j_choices == 0){
        j_choices_ret <- j_choices_mat
    } else{
        j_choices_ret <- cbind(j_choices, j_choices_mat)
    }

    if(ncol(j_choices_ret) != n_cols_j_choices + 1){
        print("ncol(j_choices_ret):")
        print(ncol(j_choices_ret))
        print("n_cols_j_choices:")
        print(n_cols_j_choices)
        stop("ncol(j_choices_ret) != n_cols_j_choices + 1")
    }

    return(list(losses_ret, j_choices_ret))
}

dataResultsLoopNoTest <- function(X_selec, y_selec, X_train, y_train, mu_train,
    R, p_max, all_cluster_feats, losses=NA, j_choices=NA, proxy_counts=NA,
    all_proxy_counts=NA){
    # Current number of columns in losses (needs to be 1 more than this later)

    if(all(is.na(losses))){
        n_cols_losses <- 0
    } else{
        n_cols_losses <- ncol(losses)
    }

    if(n_cols_losses != 0){
        if (nrow(losses) != 3*p_max){
            print("p_max:")
            print(p_max)
            print("dim(losses):")
            print(dim(losses))
            stop("nrow(losses) != 3*p_max")
        }
    }

    if(all(is.na(j_choices))){
        n_cols_j_choices <- 0
    } else{
        n_cols_j_choices <- ncol(j_choices)
    }

    if(all(is.na(proxy_counts))){
        n_cols_proxies <- 0
    } else{
        n_cols_proxies <- ncol(proxy_counts)
    }

    if(all(is.na(all_proxy_counts))){
        n_cols_all_proxies <- 0
    } else{
        n_cols_all_proxies <- ncol(all_proxy_counts)
    }

    # Possible sizes of selected sets
    j_choices_range <- 1:p_max
    # Get selected sets
    lasso_results <- getLassoSelected(X_selec, y_selec, j_choices_range)
    lasso_selected <- lasso_results[[1]]
    j_choices_range <- lasso_results[[2]]
    stab_results <- getStabSelected(X_selec, y_selec, j_choices_range, q=p_max,
        B=100)
    stab_selected <- stab_results[[1]]
    j_choices_range <- stab_results[[2]]
    gen_stab_results <- getGenStabSelected(X_selec, y_selec, j_choices_range, R=R,
        q=p_max, B=100)
    gen_stab_selected <- gen_stab_results[[1]]
    j_choices_range <- gen_stab_results[[2]]

    # Fit OLS fits for lasso selected set, stability selection selected set,
    # and generalized stability selected set
    losses_mat <- matrix(NA, nrow=3*p_max, ncol=1)
    proxies_count_mat <- matrix(NA, nrow=3*p_max, ncol=1)
    all_proxies_count_mat <- matrix(NA, nrow=3*p_max, ncol=1)
    if(!all(j_choices_range %in% 1:p_max)){
        stop("!all(j_choices_range %in% 1:p_max)")
    }
    for(j in j_choices_range){
        # losses[j, k] <- get_loss(x_train, y_train, x_test,
        #     y_test, selected_sets[[k]][[j]])
        losses_mat[j, 1] <- get_loss_no_test(X_train, y_train,
            mu_train, lasso_selected[[j]])
        # print("lasso_selected[[j]]:")
        # print(lasso_selected[[j]])
        losses_mat[j + p_max, 1] <- get_loss_no_test(X_train, y_train,
            mu_train, stab_selected[[j]])
        # print("stab_selected[[j]]:")
        # print(stab_selected[[j]])
        losses_mat[j + 2*p_max, 1] <- get_loss_no_test(X_train, y_train,
            mu_train, gen_stab_selected[[j]])
        # print("gen_stab_selected[[j]]:")
        # print(gen_stab_selected[[j]])
    }

    # Get variable numbers for proxies for real_gdp_growth
    proxy_1 <- which(colnames(data) == "real_gnp_growth")
    proxy_2 <- which(colnames(data) == "real_gdp_nfb_growth")
    proxy_nums <- c(proxy_1, proxy_2)

    for(j in j_choices_range){
        proxies_count_mat[j, 1] <- as.integer(any(proxy_nums %in% lasso_selected[[j]]))
        proxies_count_mat[j + p_max, 1] <- as.integer(any(proxy_nums %in% stab_selected[[j]]))
        proxies_count_mat[j + 2*p_max, 1] <- as.integer(any(proxy_nums %in% gen_stab_selected[[j]]))
    }

    for(j in j_choices_range){
        all_proxies_count_mat[j, 1] <- sum(lasso_selected[[j]] %in% all_cluster_feats)
        all_proxies_count_mat[j + p_max, 1] <- sum(stab_selected[[j]] %in% all_cluster_feats)
        all_proxies_count_mat[j + 2*p_max, 1] <- sum(gen_stab_selected[[j]] %in% all_cluster_feats)
    }

    # print("all_proxies_count_mat:")
    # print(all_proxies_count_mat)
    if(all(is.na(losses))){
        losses_ret <- losses_mat
    } else{
        losses_ret <- cbind(losses, losses_mat)
    }

    if(ncol(losses_ret) != n_cols_losses + 1){
        print("ncol(losses_ret):")
        print(ncol(losses_ret))
        print("n_cols_losses:")
        print(n_cols_losses)
        stop("ncol(losses_ret) != n_cols_losses + 1")
    }

    if(all(is.na(proxy_counts))){
        proxy_counts_ret <- proxies_count_mat
    } else{
        proxy_counts_ret <- cbind(proxy_counts, proxies_count_mat)
    }

    if(all(is.na(all_proxy_counts))){
        all_proxy_counts_ret <- all_proxies_count_mat
    } else{
        all_proxy_counts_ret <- cbind(all_proxy_counts, all_proxies_count_mat)
    }

    # print("all_proxy_counts_ret:")
    # print(all_proxy_counts_ret)

    if(ncol(proxy_counts_ret) != n_cols_proxies + 1){
        print("ncol(proxy_counts_ret):")
        print(ncol(proxy_counts_ret))
        print("n_cols_proxies:")
        print(n_cols_proxies)
        stop("ncol(proxy_counts_ret) != n_cols_proxies + 1")
    }

    j_choices_mat <- matrix(FALSE, nrow=3*p_max, ncol=1)
    for(i in 1:3){
        j_choices_mat[j_choices_range + p_max*(i-1), 1] <- TRUE
    }
    j_choices_mat <- as.matrix(j_choices_mat)
    if(!all(j_choices_mat[1:p_max, ] ==
        j_choices_mat[(p_max + 1):(2*p_max), ]) |
        !all(j_choices_mat[1:p_max, ] ==
        j_choices_mat[(2*p_max + 1):(3*p_max), ])){
        print(j_choices_mat)
        stop("!all.equal(j_choices_mat[1:p_max, ],
        j_choices_mat[(p_max + 1):(2*p_max), ]) |
        !all.equal(j_choices_mat[1:p_max, ],
        j_choices_mat[(2*p_max + 1):(3*p_max), ])")
    }

    if(n_cols_j_choices == 0){
        j_choices_ret <- j_choices_mat
    } else{
        j_choices_ret <- cbind(j_choices, j_choices_mat)
    }

    if(ncol(j_choices_ret) != n_cols_j_choices + 1){
        print("ncol(j_choices_ret):")
        print(ncol(j_choices_ret))
        print("n_cols_j_choices:")
        print(n_cols_j_choices)
        stop("ncol(j_choices_ret) != n_cols_j_choices + 1")
    }

    return(list(losses_ret, j_choices_ret, proxy_counts_ret,
        all_proxy_counts_ret))
}

getR <- function(snps_mat, cor_cutoff, verbose=FALSE){
    print("Calculating correlations in data...")
    t0 <- Sys.time()

    # In particular, we approximated the standardized joint distribution as 
    # multivariate Gaussian with covariance matrix estimated using the methodology 
    # of Wen and Stephens (2010), which shrinks the off- diagonal entries of the 
    # empirical covariance matrix using genetic distance information estimated from 
    # the HapMap CEU population. This approximation was used on each chromosome, 
    # and SNPs on differ- ent chromosomes were assumed to be independent. The 
    # statistic we use is the LCD. Although the data itself cannot be made 
    # available, all code is available at 
    #  http://web.stanford.edu/ ̃msesia/ software.html.

    cor_mat <- cor(snps_mat)

    zero_sds <- apply(snps_mat, 2, function(x){length(unique(x)) == 1})

    print("Number of features on this iteration with 0 SD:")

    print(sum(zero_sds))

    zero_sd_inds <- which(zero_sds)

    for(i in zero_sd_inds){
        cor_mat[, i] <- 0
        cor_mat[i, ] <- 0
        cor_mat[i, i] <- 1
    }

    print("Done! Time to calculate:")
    print(Sys.time() - t0)

    if(verbose){
        # Plot correlations

        cors <- as.data.frame(cor_mat[lower.tri(cor_mat)])

        colnames(cors) <- "Correlations"

        cor_plot <- ggplot(cors, aes(x=Correlations)) + geom_histogram()

        print(cor_plot)

        rm(cor_plot)

        print(summary(cors))

        print("Proportion of absolute correlations greater than 0.5:")
        print(sum(abs(cors$Correlations) > 0.5)/length(cors$Correlations))

        print("str(snps_mat):")
        print(str(snps_mat))
    }
    
    # print(ggcorrplot(cor_mat))

    # print(sort(abs(cor_mat[, colnames(cor_mat)=="gdp_growth"]),
    #   decreasing=TRUE))

    dist <- as.dist(1 - abs(cor_mat))
    h <- hclust(dist)
    rm(cor_mat)
    # plot(h)
    # View(data)

    # we clustered the SNPs using the estimated correlations as a similarity
    # measure with a single-linkage cutoff of 0.5, and settle for discovering 
    # important SNP clusters.

    ct <-  cutree(h, h=cor_cutoff)
    # ct
    # table(ct)

    if(verbose){
        print("Number of clusters:")

        print(length(unique(ct)))

        print("Number of SNPs:")

        print(ncol(snps_mat))

        print("Average number of SNPs per cluster:")

        print(ncol(snps_mat)/length(unique(ct)))

    }

    # Which clusters contain more than one feature?
    cluster_names <- as.integer(names(table(ct)[table(ct) > 1]))
    clustered_feats <- character()

    n_clusters <- length(cluster_names)

    if(verbose){

        for(i in 1:n_clusters){
            print(paste("Cluster ", cluster_names[i], ":", sep=""))
            clustered_feats_i <- names(ct[ct==cluster_names[i]])
            print(clustered_feats_i)
            clustered_feats <- c(clustered_feats, clustered_feats_i)
        }

        print("Unclustered features:")
        print(setdiff(colnames(data), clustered_feats))
    }

    # # Which clusters contain more than two features?
    # big_cluster_names <- as.integer(names(table(ct)[table(ct) > 2]))

    # Create R matrix
    p <- ncol(snps_mat)
    R <- diag(1, p)
    colnames(R) <- colnames(snps_mat)
    rownames(R) <- colnames(snps_mat)

    for(i in 1:n_clusters){
        cluster_i_vars <- names(ct[ct==cluster_names[i]])
        if(length(cluster_i_vars) <= 1){
            stop("length(cluster_i_vars) <= 1")
        }
        for(j in 1:(length(cluster_i_vars) - 1)){
            for(k in j:length(cluster_i_vars)){
                R[cluster_i_vars[j], cluster_i_vars[k]] <- 1
                R[cluster_i_vars[k], cluster_i_vars[j]] <- 1
            }
        }
    }

    if(!identical(R, t(R))){
        stop("R is not symmetric")
    }
    if(any(diag(R) != 1)){
        stop("R does not contain all entries on the diagonal equal to 1.")
    }
    return(R)
}



getRankedFeatures <- function(df, method, j_choices, var_names, response_name, 
    R=NA){

    t2 <- Sys.time()

    # List of selected sets to return
    selected_sets <- list()
    if(length(j_choices) == 0){
        warning("j_choices is empty")
        return(list(selected_sets, j_choices))
    }

    if(!(response_name %in% colnames(df))){
        stop("!(response_name %in% colnames(df))")
    }

    # j's to elminate from j_choices
    j_choices_to_eliminate <- integer()

    x <- as.matrix(df[, colnames(df) != response_name])
    y <- df[, response_name]

    if(method=="lasso"){
        fit <- glmnet(x=x, y=y, family="gaussian", alpha=1, nlambda=2000,
            lambda.min.ratio=0.1)
        selected <- unique(predict(fit, type="nonzero"))

        lasso_sets_results <- extractLassoSets(selected, j_choices,
            j_choices_to_eliminate, var_names)

        selected_sets <- lasso_sets_results$selected_sets
        j_choices_to_eliminate <- lasso_sets_results$j_choices_to_eliminate

        rm(lasso_sets_results)
        
    } else if(method=="lasso_proto"){

        if(any(is.na(R))){
            stop(paste("must provide R for method", method))
        }

        protolasso_fit <- protolasso(x, y, R, var_names=var_names)

        protolasso_results <- extractProtolassoSets(protolasso_fit$selected_sets,
            j_choices, j_choices_to_eliminate, var_names=var_names)

        selected_sets <- protolasso_results$selected_sets
        j_choices_to_eliminate <- protolasso_results$j_choices_to_eliminate

        rm(protolasso_fit)
        rm(protolasso_results)
        
    } else if(method=="SS_SS"){
        # Stability selection
        B <- 100
        n <- nrow(x)

        stopifnot(length(y) == n)
        # Determine lambda: do lasso with cross-validation on full data sample.
        # Cross-validated model size will be the size we use for stability
        # selection.
        inds_size <- sample(1:n, floor(n/2))

        size_results <- cv.glmnet(x=x[inds_size, ], y=y[inds_size],
            parallel=TRUE, family="gaussian")

        fit <- css(X=x, y=y, lambda=size_results$lambda.min,
            B=B, sampling_type="SS")

        GSS_sets_results <- extractGSSSets(fit$selected,
            getTiedWithNext(names(fit$selected_clusts), fit$clus_sel_props),
            j_choices, j_choices_to_eliminate, var_names=var_names)

        selected_sets <- GSS_sets_results$selected_sets
        j_choices_to_eliminate <- GSS_sets_results$j_choices_to_eliminate

        rm(GSS_sets_results)

    } else if(method=="SS_GSS"){
        # Sparse generalized stability selection (Shah and Samworth version)
        B <- 50
        n <- nrow(x)

        if(any(is.na(R))){
            stop(paste("must provide R for method", method))
        }

        stopifnot(length(y) == n)
        # Determine lambda: do lasso with cross-validation on full data sample.
        # Cross-validated model size will be the size we use for stability
        # selection.
        inds_size <- sample(1:n, floor(n/2))

        clusters <- getClustersFromR(R)

        size_results <- cv.glmnet(x=x[inds_size, ], y=y[inds_size],
            parallel=TRUE, family="gaussian")

        fit <- css(X=x, y=y, lambda=size_results$lambda.min, clusters=clusters,
            B=B, sampling_type="SS", weighing="sparse")

        GSS_sets_results <- extractGSSSets(fit$selected, getTiedWithNext(names(fit$selected_clusts), fit$clus_sel_props),
            j_choices, j_choices_to_eliminate, var_names=var_names)

        selected_sets <- GSS_sets_results$selected_sets
        j_choices_to_eliminate <- GSS_sets_results$j_choices_to_eliminate

        rm(GSS_sets_results)

    } else if(method=="SS_GSS_avg"){
        # Weighted averaged generalized stability selection (Shah and Samworth
        # version)
        B <- 100

        if(any(is.na(R))){
            stop(paste("must provide R for method", method))
        }

        n <- nrow(x)
        if(length(y) != n){
            stop("length(y) != n")
        }
        # Determine q: do lasso with cross-validation on full data sample.
        # Cross-validated model size will be the size we use for stability
        # selection.
        inds_size <- sample(1:n, floor(n/2))
        size_results <- cv.glmnet(x=x[inds_size, ], y=y[inds_size],
            parallel=TRUE, family="gaussian")

        clusters <- getClustersFromR(R)

        fit <- css(X=x, y=y, lambda=size_results$lambda.min, clusters=clusters,
            B=B, sampling_type="SS", weighting="weighted_avg")


        GSS_avg_results <- extractGSS_avgSets(fit$selected, fit$to_avg,
            fit$avg_feats, fit$weights, getTiedWithNext(names(fit$selected_clusts), fit$clus_sel_props), j_choices,
            j_choices_to_eliminate, var_names=var_names)

        selected_sets <- GSS_avg_results$selected_sets
        j_choices_to_eliminate <- GSS_avg_results$j_choices_to_eliminate
        to_avg_list <- GSS_avg_results$to_avg_list
        avg_feats_list <- GSS_avg_results$avg_feats_list
        weights_list <- GSS_avg_results$weights_list

        rm(GSS_avg_results)

    } else if(method=="SS_GSS_avg_unwt"){
        # Simple averaged generalized stability selection (Shah and Samworth
        # version)
        # cutoff <- 0.52
        B <- 100

        n <- nrow(x)

        if(length(y) != n){
            stop("length(y) != n")
        }

        if(any(is.na(R))){
            stop(paste("must provide R for method", method))
        }

        # Determine lambda: do lasso with cross-validation on full data sample.
        # Cross-validated model size will be the size we use for stability
        # selection.
        inds_size <- sample(1:n, floor(n/2))
        size_results <- cv.glmnet(x=x[inds_size, ], y=y[inds_size],
            parallel=TRUE, family="gaussian")

        clusters <- getClustersFromR(R)

        fit <- css(X=x, y=y, lambda=size_results$lambda.min,
            clusters=clusters, B=B, sampling_type="SS", weighting="simple_avg")

        GSS_avg_results <- extractGSS_avgSets(fit$selected, fit$to_avg,
            fit$avg_feats, fit$weights,
            getTiedWithNext(names(fit$selected_clusts), fit$clus_sel_props),
            j_choices,
            j_choices_to_eliminate, var_names=var_names)

        selected_sets <- GSS_avg_results$selected_sets
        j_choices_to_eliminate <- GSS_avg_results$j_choices_to_eliminate
        to_avg_list <- GSS_avg_results$to_avg_list
        avg_feats_list <- GSS_avg_results$avg_feats_list
        weights_list <- GSS_avg_results$weights_list

        rm(GSS_avg_results)

    } else if(method=="BRVZ_avg_unwt"){

        if(any(is.na(R))){
            stop(paste("must provide R for method", method))
        }

        BRVZ_fit <- clusterRepLasso(x, y, R, var_names=var_names)

        # selected <- unique(BRVZ_results$selected)
        # j_choices_to_eliminate <- BRVZ_results$j_choices_to_eliminate

        BRVZ_results <- extractBRVZSets(BRVZ_fit$selected_sets,
            BRVZ_fit$to_avg_list, BRVZ_fit$avg_feats_list,
            BRVZ_fit$weights_list, j_choices, j_choices_to_eliminate,
            var_names=var_names)

        selected_sets <- BRVZ_results$selected_sets
        j_choices_to_eliminate <- BRVZ_results$j_choices_to_eliminate
        to_avg_list <- BRVZ_results$to_avg_list
        avg_feats_list <- BRVZ_results$avg_feats_list
        weights_list <- BRVZ_results$weights_list

        rm(BRVZ_results)
        rm(BRVZ_fit)

    } else{
        message <- paste("No feature extraction method exists for method",
            nameMap(method), ".")
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

    if(length(selected_sets) < max(j_choices)){
        print("method:")
        print(method)
        print("length(selected_sets):")
        print(length(selected_sets))
        print("max(j_choices):")
        print(max(j_choices))
        stop("length(selected_sets) < max(j_choices)")
    }

    if(!(method %in% c("SS_GSS_avg", "SS_GSS_avg_unwt", "BRVZ_avg_unwt"))){
        return(list(selected_sets=selected_sets, j_choices=j_choices))
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
    # will contain a character vector of the names of the features from the
    # cluster of which feature k is a member.
    #
    # weights_list: A list whose elements are lists of length j. If
    # feature k in 
    # 1:j is in a cluster, the kth entry of the jth vector of
    # weights_list will
    # be the weights to use (in the same order as the jth entry of
    # avg_feats_list).

    return(list(selected_sets=selected_sets, j_choices=j_choices,
        to_avg_list=to_avg_list, avg_feats_list=avg_feats_list,
        weights_list=weights_list))

}

getLossPlant <- function(df_train, df_test, selected, var_names, response_name,
    average=FALSE, to_avg=NA, avg_feats=NA, weights=NA) {
    # Inputs:
    #
    # selected: a character vector of length j containing the names of the
    # selected features
    #
    # average: if TRUE, will do an averaging of cluster members.
    #
    # to_avg: must be provided if average=TRUE. A logical vector of the same
    # length as selected. If feature j is in a cluster, the jth entry of 
    # to_avg will be TRUE. 
    #
    # avg_feats: must be provided if average=TRUE. A list of the same length as
    # selected . If feature j is in a
    # cluster, the jth entry of avg_feats will contain a namd integer vector
    # of the indices of the features from the cluster of which feature j is a
    # member.
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

    if(!(response_name %in% colnames(df_train))){
        stop("!(response_name %in% colnames(df_train))")
    }

    if(!(response_name %in% colnames(df_test))){
        stop("!(response_name %in% colnames(df_test))")
    }

    n_train <- nrow(df_train)
    n_test <- nrow(df_test)

    x_train <- df_train[, colnames(df_train) != response_name]
    if(ncol(x_train) != length(var_names)){
        stop("ncol(x_train) != length(var_names)")
    }
    colnames(x_train) <- var_names
    y_train <- df_train[, response_name]

    x_test <- df_test[, colnames(df_test) != response_name]
    if(ncol(x_test) != length(var_names)){
        stop("ncol(x_test) != length(var_names)")
    }
    colnames(x_test) <- var_names
    y_test <- df_test[, response_name]

    p <- ncol(x_train)

    if(p != ncol(x_test)){
        stop("p != ncol(x_test)")
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
            stop("length(to_avg) > ncol(x_train)")
        }
        if(!all.equal(lengths(avg_feats), lengths(weights))){
            stop("!all.equal(lengths(avg_feats), lengths(weights))")
        }
    }

    if(!all(selected %in% var_names)){
        stop("!all(selected %in% var_names)")
    }

    data_train <- createDesignForLoss2(x_train, y_train, selected, var_names,
        average, to_avg, avg_feats, weights)

    colnames(data_train)[2:ncol(data_train)] <- selected

    data_test <- createDesignForLoss2(x_test, y_test, selected, var_names,
        average, to_avg, avg_feats, weights)

    colnames(data_test)[2:ncol(data_test)] <- selected

    if(any(colnames(data_test)[2:ncol(data_test)] !=
        colnames(data_train)[2:ncol(data_train)])){
        colnames(data_test)[2:ncol(data_test)] <-
            colnames(data_train)[2:ncol(data_train)]
    }

    model <- suppressWarnings(lm(y~., data_train))
    loss <- Metrics::mse(actual=y_test,
        predicted=suppressWarnings(predict(model, newdata=data_test)))
    return(loss)
}

createDesignForLoss2 <- function(X, y, selected, var_names, average=FALSE,
    to_avg=NA, avg_feats=NA, weights=NA){
    # Inputs:
    #
    # selected: a character vector of length j containing the names of the
    # selected features
    # 
    # average: if TRUE, will do an averaging of cluster members.
    #
    # to_avg: must be provided if average=TRUE. A logical vector of the same
    # length as selected. If feature j is in a cluster, the jth entry of 
    # to_avg will be TRUE. 
    #
    # avg_feats: must be provided if average=TRUE. A list of length at most the
    # same length as selected . If feature j is in a
    # cluster, the jth entry of avg_feats will contain a character vector
    # containing the names of the features from the cluster of which feature j 
    # is a member.
    #
    # weights: must be provided if average=TRUE. A list the same length as
    # avg_feats. If feature j is in a cluster, the jth entry of weights will
    # be the weights to use (in the same order as the jth entry of
    # avg_feats).

    if(average){
        if(all(is.na(to_avg))){
            stop("If average==TRUE, must provide to_avg to get_loss")
        }
        if(any(to_avg)){
            if(all(is.na(avg_feats)) | all(is.na(weights))){
                stop("If average==TRUE, must provide avg_feats and weights to get_loss")
            }
        }
        if(length(to_avg) > ncol(X)){
            stop("length(to_avg) > ncol(X)")
        }
        if(!all.equal(lengths(avg_feats), lengths(weights))){
            stop("!all.equal(lengths(avg_feats), lengths(weights))")
        }
    }

    p_selected <- length(selected)

    if(p_selected > ncol(X)){
        stop("p_selected > ncol(X)")
    }

    if(!all(selected %in% var_names)){
        stop("!all(selected %in% var_names)")
    }

    if(!average){
        X_df <- X[, var_names %in% selected]
    } else{
        X_df <- matrix(0, nrow(X), p_selected)
        colnames(X_df) <- character(p_selected)
        for(j in 1:p_selected){
            if(to_avg[j]){
                # Features to average
                weights_j <- weights[[j]]
                if(abs(sum(weights_j) - 1) > 10^(-6)){
                    stop("abs(sum(weights_j) - 1) > 10^(-6)")
                }
                if(length(avg_feats[[j]]) != length(weights_j)){
                    stop("length(avg_feats[[j]]) != length(weights_j)")
                }
                # print("avg_feats[[j]]:")
                # print(avg_feats[[j]])
                # print("weights_j:")
                # print(weights_j)
                # print("str(X):")
                # print(str(X))
                # print("str(X[, avg_feats[[j]]]):")
                # print(str(X[, avg_feats[[j]]]))
                X_df[, j] <- as.matrix(X[, var_names %in% avg_feats[[j]]]) %*%
                    weights_j
                if(p_selected > 1){
                    # print("str(X_df):")
                    # print(str(X_df))
                    # print(paste("X", selected[j], sep=""))
                    # print("colnames(X_df):")
                    # print(colnames(X_df))
                    # print("j:")
                    # print(j)
                    colnames(X_df)[j] <- selected[j]
                }
            } else{
                X_df[, j] <- X[, var_names == selected[j]]
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
        colnames(df)[2] <- selected
    }

    return(df)
}

getErrorDf3 <- function(losses, names, j_choices_mat, p_max_plots, J,
    coarseness){

    if(any(losses < 0)){
        stop("any(losses < 0)")
    }

    n_methods <- length(names)
    n_draws <- ncol(losses)

    if(is.na(J)){
        stop("is.na(J)")
    }

    if(n_draws > 1){
        p_max <- nrow(losses)/n_methods

        if(p_max != round(p_max)){
            stop("p_max != round(p_max)")
        }

        n_losses <- nrow(losses)/coarseness

        if(n_losses != round(n_losses)){
            stop("n_losses != round(n_losses)")
        }

        any_losses <- logical(n_losses)

        losses_vec <- as.numeric(rep(NA, n_losses))

        for(i in 1:n_losses){
            inds_i <- (coarseness*(i - 1) + 1):(coarseness*i)
            any_losses[i] <- sum(losses[inds_i, ]) > 0
            # Confirm this is working the way I'm expecting
            if(any_losses[i] != any(j_choices_mat[inds_i, ])){
                stop("any_losses[i] != any(j_choices_mat[inds_i, ])")
            }
        }





        # # j_choices_logical <- rowSums(j_choices_mat[(J*(n_methods - 1) + 1):
        # #     (J*n_methods), ]) > 0
        # j_choices_logical <- rowSums(j_choices_mat) > 0
        # # j_choices <- which(j_choices_logical)
        # # p_max <- max(j_choices)
        # # print("j_choices:")
        # # print(j_choices)
        # # row_inds <- which(rep(j_choices_logical, n_methods))
        # row_inds <- which(j_choices_logical)
        
        # # n_losses <- ceiling(sum(j_choices_logical)/coarseness)*n_methods
        # n_losses <- length(row_inds)
        # # print("sum(j_choices_logical):")
        # # print(sum(j_choices_logical))
        # # print("n_methods:")
        # # print(n_methods)
        # # print("coarseness:")
        # # print(coarseness)
        # # print("n_losses:")
        # # print(n_losses)
        # losses_vec <- numeric(n_losses)
        
        # # if(ceiling(length(j_choices)/coarseness)*n_methods != n_losses){
        # #     print("length(j_choices):")
        # #     print(length(j_choices))
        # #     print("n_methods:")
        # #     print(n_methods)
        # #     print("coarseness:")
        # #     print(coarseness)
        # #     stop("ceiling(length(j_choices)/coarseness)*n_methods != n_losses")
        # # }
        # # if(length(row_inds) != n_losses){
        # #     stop("length(row_inds) != n_losses")
        # # }


        for(i in 1:n_losses){
            if(any_losses[i]){
                row_inds_i <- (coarseness*(i - 1) + 1):(coarseness*i)
                # print("row_inds_i:")
                # print(row_inds_i)
                cols_i <- j_choices_mat[row_inds_i, ]
                if(!identical(dim(losses[row_inds_i, ]), dim(cols_i))){
                    stop("!identical(dim(losses[row_inds_i, ]), dim(cols_i))")
                }
                # if(coarseness > 1){
                #     for(j in 2:coarseness){
                #         cols_i <- c(cols_i, j_choices_mat[row_inds_i[j], ])
                #     }
                # }
                # print("cols_i:")
                # print(cols_i)
                losses_vec[i] <- mean(losses[row_inds_i, ][cols_i])
            }
            
        }
    } else{
        stop("haven't revised code for this section yet")
        row_inds <- which(j_choices_mat)
        n_losses <- ceiling(sum(j_choices_mat)/coarseness)
        losses_vec <- numeric(n_losses)
        j_choices <- which(j_choices_mat[1:J, 1])
        if(ceiling(length(j_choices)*n_methods/coarseness) != n_losses){
            stop("ceiling(length(j_choices)*n_methods/coarseness) != n_losses")
        }
        # if(length(row_inds) != n_losses){
        #     stop("length(row_inds) != n_losses")
        # }

        for(i in 1:n_losses){
            row_inds_i <- row_inds[((i-1)*coarseness + 1):min(i*coarseness,
                nrow(row_inds))]
            losses_vec[i] <- mean(losses[row_inds_i, 1])
        }

        # losses_vec <- losses[row_inds, 1]
    }

    if(n_losses/n_methods != round(n_losses/n_methods)){
        print("n_losses:")
        print(n_losses)
        print("n_methods:")
        print(n_methods)
        stop("n_losses/n_methods != round(n_losses/n_methods)")
    }

    # j_choices_min <- coarseness*(0:(n_losses - 1)) + 1
    # j_choices_min <- j_choices_min %% J
    # j_choices_max <- coarseness*(1:n_losses)
    # j_choices_max <- j_choices_max %% J

    j_choices_min <- coarseness*(0:(J/coarseness - 1)) + 1
    # j_choices_min <- j_choices_min %% J
    j_choices_max <- coarseness*(1:(J/coarseness))
    # j_choices_max <- j_choices_max %% J

    # j_choices_min <- integer(n_losses/n_methods)
    # j_choices_max <- integer(n_losses/n_methods)
    
    # for(i in 1:(n_losses/n_methods)){
    #     row_inds_i <- row_inds[((i-1)*coarseness + 1):min(i*coarseness,
    #         nrow(row_inds))]
    #     j_choices_min[i] <- min(row_inds_i)
    #     j_choices_max[i] <- max(row_inds_i)
    # }

    j_choices_display <- (j_choices_min + j_choices_max)/2

    if(length(j_choices_display)*n_methods != length(losses_vec)){
        print("length(j_choices_display):")
        print(length(j_choices_display))
        print("n_methods:")
        print(J/coarseness)
        print("length(losses_vec)):")
        print(length(losses_vec))
        stop("length(j_choices_display)*n_methods != length(losses_vec)")
    }

    if(length(j_choices_display) != length(unique(j_choices_display))){
        stop("length(j_choices_display) != length(unique(j_choices_display))")
    }

    if(length(losses_vec) != n_methods*J/coarseness){
        print("length(losses_vec):")
        print(length(losses_vec))
        print("n_methods:")
        print(n_methods)
        print("J:")
        print(J)
        stop("length(losses_vec) != n_methods*J/coarseness")
    }

    if(J/coarseness != round(J/coarseness)){
        stop("J/coarseness != round(J/coarseness)")
    }

    # row_inds <- j_choices %% J
    # row_inds[row_inds == 0] <- J

    df_gg <- data.frame(rep(j_choices_display, times=n_methods), losses_vec,
        rep(names, each=J/coarseness))

    colnames(df_gg) <- c("Model_Size", "MSE", "Method")

    if(any(duplicated(df_gg[, c("Model_Size", "Method")]))){
        print("df_gg[duplicated(df_gg[, c(Model_Size, Method)]), ]:")
        print(df_gg[duplicated(df_gg[, c("Model_Size", "Method")]), ])
        stop("any(duplicated(df_gg[, c(Model_Size, Method)]))")
    }

    df_gg <- df_gg[complete.cases(df_gg), ]

    if(!is.na(p_max_plots)){
        df_gg <- df_gg[df_gg$Model_Size <= p_max_plots, ]
    }

    if(any(duplicated(df_gg[, c("Model_Size", "Method")]))){
        print("df_gg[duplicated(df_gg[, c(Model_Size, Method)]), ]:")
        print(df_gg[duplicated(df_gg[, c(Model_Size, Method)]), ])
        stop("any(duplicated(df_gg[, c(Model_Size, Method)]))")
    }
    

    # print("head(df_gg):")
    # print(head(df_gg))

    # print("summary(df_gg):")
    # print(summary(df_gg))

    # print("str(df_gg):")
    # print(str(df_gg))

    # print("unique(df_gg$Method):")
    # print(unique(df_gg$Method))

    # Make labels display-friendly
    df_gg <- makeDfMethodsDisplay3(df_gg)

    return(df_gg)
}

createLossesPlot3 <- function(losses, names, j_choices_mat, p_max_plots=NA,
    J=NA, legend=TRUE, coarseness=1){
    n_methods <- length(names)
    n_draws <- ncol(losses)

    # If J is not provided, figure out on own
    if(length(J) != 1){
        stop("length(J) != 1")
    }
    if(is.na(J)){
        J <- nrow(losses)/n_methods
    }
    if(J != round(J)){
        print(nrow(losses))
        print(n_methods)
        print(J)
        stop("J != round(J)")
    }
    if(!is.na(p_max_plots)){
        if(p_max_plots/coarseness != round(p_max_plots/coarseness)){
            stop("p_max_plots/coarseness != round(p_max_plots/coarseness)")
        }

    } else{
        if(J/coarseness != round(J/coarseness)){
            stop("J/coarseness != round(J/coarseness)")
        }
    }
    if(any(losses < 0)){
        stop("any(losses < 0)")
    }
    if(n_draws != ncol(j_choices_mat) | nrow(losses) != nrow(j_choices_mat)){
        print(dim(losses))
        print(dim(j_choices_mat))
        stop("ncol(losses) != ncol(j_choices_mat) | nrow(losses) != nrow(j_choices_mat)")
    }

    if(length(p_max_plots) != 1){
        stop("length(p_max_plots) != 1")
    }

    df_gg <- getErrorDf3(losses, names, j_choices_mat, p_max_plots, J,
        coarseness)

    # print("head(df_gg):")
    # print(head(df_gg))

    # print("summary(df_gg):")
    # print(summary(df_gg))

    # print("str(df_gg):")
    # print(str(df_gg))

    # print("unique(df_gg$Method):")
    # print(unique(df_gg$Method))

    # Get max ranking value for integer labels on horizontal axis
    max_rank <- max(df_gg$Model_Size)
    if((max_rank %% 2) != 0){
        max_rank <- max_rank + 1
    }

    # print("max_rank:")
    # print(max_rank)

    plot <- ggplot(df_gg, aes(x=Model_Size, y=MSE, color=Method,
        shape=Method)) + scale_shape_manual(values=1:n_methods) + 
        # geom_point() +
        geom_point(size=2.5, alpha=1) +
        xlab("No. Fitted Coefficients")+
        ylim(0.03, .05)
        # + scale_x_continuous(breaks=seq(2, max_rank, by=2))

    if(!legend){
        plot <- plot + theme(legend.position="none")
    }

    return(plot)
}

makeDfMethodsDisplay3 <- function(df){
    # Takes in data.frame with computer-friendly labels and returns 
    # data.frame with display-friendly labels
    if(!("Method" %in% colnames(df))){
        stop("!(Method %in% colnames(df))")
    }
    if(ncol(df) != 3){
        stop("ncol(df) != 3 (before modifications)")
    }
    labels <- as.factor(nameMap(as.character(df$Method)))
    labels <- droplevels(labels)
    df$Method <- labels
    if(!all(levels(df$Method) %in% c("Lasso", "Stability Selection",
        "Sparse CSS", "Weighted Averaged CSS",
        "Simple Averaged CSS", "Cluster Representative Lasso",
        "Protolasso"))){
        stop("!all(levels(df$Method) %in% c(Lasso, Stability Selection,
            Sparse CSS, Weighted Averaged CSS,
            Simple Averaged CSS, Cluster Representative Lasso, Protolasso
))")
    }
    if(ncol(df) != 3){
        stop("ncol(df) != 3 (after modifications)")
    }
    return(df)
}

# nameMap3 <- function(sys_name){
#     # Takes in computer-friendly name for method and outputs display name
#     if(!is.character(sys_name)){
#         stop("!is.character(sys_name)")
#     }
#     ret <- sys_name
#     ret[sys_name %in%  c("lasso", "lasso_random")] <- "Lasso"
#     ret[sys_name %in%  c("SS_SS")] <- "Stability Selection"
#     ret[sys_name %in%  c("SS_GSS")] <- "Sparse CSS"
#     ret[sys_name %in%  c("SS_GSS_avg")] <- "Weighted Averaged CSS"
#     ret[sys_name %in%  c("SS_GSS_avg_unwt")] <- "Simple Averaged CSS"
#     return(ret)
# }

getNSBDf3 <- function(stab_mets, coarseness=1, plot_errors=FALSE, lowers=NA,
    uppers=NA){

    n_methods <- ncol(stab_mets)
    p_max <- nrow(stab_mets)*coarseness
    j_choices_max <- coarseness*(1:nrow(stab_mets))
    j_choices_min <- 1 + coarseness*(0:(nrow(stab_mets) - 1))
    j_choices <- (j_choices_max + j_choices_min)/2
    if(length(j_choices) != nrow(stab_mets)){
        stop("length(j_choices) != nrow(stab_mets)")
    }
    if(any(j_choices != round(j_choices))){
        stop("any(j_choices != round(j_choices))")
    }
    # j_choices <- 1:p_max
    names <- colnames(stab_mets)

    # # Eliminate any rows with all NA entries
    # all_na_rows <- apply(stab_mets, 1, function(x){all(is.na(x))})
    # j_choices_to_eliminate <- j_choices[all_na_rows]
    # j_choices <- setdiff(j_choices, j_choices_to_eliminate)
    # stab_mets <- stab_mets[!all_na_rows, ]
    # if(plot_errors){
    #     lowers <- lowers[!all_na_rows, ]
    #     uppers <- uppers[!all_na_rows, ]
    # }

    # # Eliminate any rows with any NA entries
    # any_na_rows <- apply(stab_mets, 1, function(x){any(is.na(x))})
    # j_choices_to_eliminate <- j_choices[any_na_rows]
    # j_choices <- setdiff(j_choices, j_choices_to_eliminate)
    # stab_mets <- stab_mets[!any_na_rows, ]
    # if(plot_errors){
    #     lowers <- lowers[!any_na_rows, ]
    #     uppers <- uppers[!any_na_rows, ]
    # }

    # # Now get j_choices_display

    # if(coarseness > 1){

    # } else{
    #     j_choices_display <- j_choices
    # }

    df_gg <- data.frame(rep(j_choices, n_methods), as.vector(stab_mets),
        rep(names, each=length(j_choices)))

    colnames(df_gg) <- c("Model_Size", "NSBStability", "Method")

    df_gg <- df_gg[complete.cases(df_gg), ]

    if(any(duplicated(df_gg[, c("Model_Size", "Method")]))){
        print("df_gg[duplicated(df_gg[, c(Model_Size, Method)]), ]:")
        print(df_gg[duplicated(df_gg[, c(Model_Size, Method)]), ])
        stop("any(duplicated(df_gg[, c(Model_Size, Method)]))")
    }

    # Make labels display-friendly
    df_gg <- makeDfMethodsDisplay3(df_gg)

    return(list(df_gg=df_gg, lowers=lowers, uppers=uppers))
}

createNSBStabPlot3 <- function(stab_mets, coarseness=1, legend=TRUE,
    plot_errors=FALSE, lowers=NA, uppers=NA){
    require(ggplot2)

    n_methods <- ncol(stab_mets)
    p_max <- nrow(stab_mets)*coarseness

    NSBDf3_res <- getNSBDf3(stab_mets=stab_mets,
        coarseness=coarseness, plot_errors=plot_errors, lowers=lowers,
        uppers=uppers)

    df_gg <- NSBDf3_res$df_gg
    lowers <- NSBDf3_res$lowers
    uppers <- NSBDf3_res$uppers

    rm(NSBDf3_res)

    # Get max ranking value for integer labels on horizontal axis
    max_rank <- max(df_gg$Model_Size)
    if((max_rank %% 2) != 0){
        max_rank <- max_rank + 1
    }

    if(plot_errors){
        # Add in error bounds
        lower_bounds <- as.vector(lowers)
        upper_bounds <- as.vector(uppers)
        if(length(lower_bounds) != nrow(df_gg) |
            length(upper_bounds) != nrow(df_gg)){
            print("length(upper_bounds):")
            print(length(upper_bounds))
            print("length(upper_bounds):")
            print(length(upper_bounds))
            print("nrow(df_gg):")
            print(nrow(df_gg))
            stop("length(lower_bounds) != nrow(df_gg) | length(upper_bounds) != nrow(df_gg)")
        }
        df_gg <- data.frame(df_gg, lower_bounds, upper_bounds)
        colnames(df_gg)[(ncol(df_gg) - 1):ncol(df_gg)] <- c("Lower", "Upper")
    }

    # # Only include lasso, stability selection, and generalized stability 
    # # selection

    # df_gg <- df_gg[df_gg$Method %in% c("Lasso", "Stability Selection",
    #     "Sparse CSS"), ]

    plot <- ggplot(df_gg, aes(x=Model_Size, y=NSBStability, color=Method,
        shape=Method)) + scale_shape_manual(values=1:n_methods) +
        geom_point(size=2.5, alpha=1) + xlab("No. Fitted Coefficients") +
        ylab("NSB Stability") + scale_x_continuous(breaks=seq(floor(p_max/10),
            max_rank, by=floor(p_max/10)))

    if(plot_errors){
        plot <- plot + geom_errorbar(aes(ymin = Lower, ymax = Upper),
            width = 0.5)
    }

    if(!legend){
        plot <- plot + theme(legend.position="none")
    }

    return(plot)
}

createStabMSEPlot3 <- function(losses, stab_mets, names, j_choices_mat,
    p_max_plots=NA, J=NA, legend=TRUE, coarseness=1, plot_errors=FALSE,
    lowers=NA, uppers=NA, names_to_omit=NA){

    require(dplyr)

    n_methods <- ncol(stab_mets)
    p_max <- nrow(stab_mets)*coarseness

    if(length(names) != n_methods){
        stop("length(names) != ncol(stab_mets)")
    }
    n_draws <- ncol(losses)

     # If J is not provided, figure out on own
    if(length(J) > 1){
        stop("length(J) > 1")
    }
    if(is.na(J)){
        J <- nrow(losses)/n_methods
    }
    if(J != round(J)){
        print(nrow(losses))
        print(n_methods)
        print(J)
        stop("J != round(J)")
    }
    if(!is.na(p_max_plots)){
        if(p_max_plots/coarseness != round(p_max_plots/coarseness)){
            stop("p_max_plots/coarseness != round(p_max_plots/coarseness)")
        }

    } else{
        if(J/coarseness != round(J/coarseness)){
            stop("J/coarseness != round(J/coarseness)")
        }
    }
    if(any(losses < 0)){
        stop("any(losses < 0)")
    }
    if(n_draws != ncol(j_choices_mat) | nrow(losses) != nrow(j_choices_mat)){
        print(dim(losses))
        print(dim(j_choices_mat))
        stop("ncol(losses) != ncol(j_choices_mat) | nrow(losses) != nrow(j_choices_mat)")
    }

    if(length(p_max_plots) != 1){
        stop("length(p_max_plots) != 1")
    }

    NSBDf3_res <- getNSBDf3(stab_mets=stab_mets,
        coarseness=coarseness, plot_errors=plot_errors, lowers=lowers,
        uppers=uppers)

    df_gg_stab <- NSBDf3_res$df_gg
    lowers <- NSBDf3_res$lowers
    uppers <- NSBDf3_res$uppers

    # print("head(df_gg_stab):")
    # print(head(df_gg_stab))

    # print(df_gg_stab[df_gg_stab$Method=="Simple Averaged CSS", ])

    rm(NSBDf3_res)

    df_gg_mse <- getErrorDf3(losses, names, j_choices_mat, p_max_plots, J,
        coarseness)

    # print("head(df_gg_mse):")
    # print(head(df_gg_mse))

    # print(df_gg_mse[df_gg_mse$Method=="Simple Averaged CSS", ])

    df_gg <- left_join(df_gg_stab, df_gg_mse, by=c("Model_Size", "Method"))

    if(any(!is.na(names_to_omit))){
        df_gg <- df_gg[!(df_gg$Method %in% names_to_omit), ]
        df_gg <- df_gg[!(df_gg$Method %in% nameMap(names_to_omit)), ]
    }

    # print("head(df_gg):")
    # print(head(df_gg))
    # print(df_gg[df_gg$Method=="Simple Averaged CSS", ])

    plot <- ggplot(df_gg, aes(x=NSBStability, y=MSE,
        color=Method, shape=Method)) + scale_shape_manual(values=1:n_methods) +
        suppressWarnings(geom_point(size=2.5, alpha=1)) + ylim(0.03, .05) +
        xlab("NSB Stability")

    if(!legend){
        plot <- plot + theme(legend.position="none")
    }

    return(plot)

}




print("simFunctions.R loaded correctly")