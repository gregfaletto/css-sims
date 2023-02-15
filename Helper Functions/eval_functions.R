## @knitr metrics

# Screening Criterion 4 (for linear models)

# This is the screening criterion.
sc4 <- function(selected, X.mat, proj_x.mu){
  # This calculates || proj_col(X_A) mu ||^2 / ||proj_col(X) mu ||^2
  #
     # # Get digest for this feature set
     # digest <- digest(selected)
     # # Check to see if we've already calculated the screening criterion for
     # # this feature set; if so, return the calculated value.
     # if(digest %in% digest.matrix[, 1]){
     #      index <- digest.matrix[, 1] == digest
     #      return(digest.matrix[index, 2])
     # }
     # # If not, calculate and then add to digest.matrix
     
     if(length(selected)==0){
          return(0)
     }
     #  # get the projection of mu onto the column space of X_A
     #   proj_xA.mu <- X.mat[, selected] %*% ginv(X.mat[, selected])%*% proj_x.mu
     #   # sv <- svd(X.mat[, selected])

     ############## Instead of inverting matrix, use linear regression package
     if(length(selected)==1){
          beta_hat_fitted <- lm.fit(as.matrix(X.mat[, selected]), proj_x.mu,
               tol=1e-12)$coefficients
     } else{
          beta_hat_fitted <- lm.fit(X.mat[, selected], proj_x.mu,
               tol=1e-12)$coefficients
     }

     # Check for possible problems
     if(!is.numeric(beta_hat_fitted)){
          print("beta_hat_fitted is not a vector.")
          print(beta_hat_fitted)
     }
     if(length(beta_hat_fitted) != length(selected)){
          print("length mismatch")
          print(length(beta_hat_fitted))
          print(length(selected))
          print(beta_hat_fitted)
     # Calculate the projection of mu onto the column space of X_A
     }
     if(length(selected)==1){
          proj_xA.mu <- X.mat[, selected]*beta_hat_fitted
     }else{
          proj_xA.mu <- X.mat[, selected] %*% beta_hat_fitted
     }
     
     # find sum of squares
     # proj_xA.mu.squared <- t(proj_xA.mu) %*% proj_xA.mu
     proj_xA.mu.squared <- sum(proj_xA.mu^2)
     # Because the true model is linear, the linear projection of
     # the true mu onto the column space of X is simply X*beta.)
     proj_x.mu.squared <- sum(proj_x.mu^2) #t(mu)%*% (mu)
    
     # Screening Criterion 4 (for linear models)
     sc4 <- proj_xA.mu.squared/proj_x.mu.squared

     # # Add to digest table
     # digest.matrix <- rbind(digest.matrix, c(digest, sc4))
     # Return screening criterion and digest table
     return(sc4)
}

sc4linear <- new_metric("sc4linear", "Screening Criterion 4 for linear models",
     metric = function(model, out) {
       # as written, this assumes that mu and proj_x.mu are the same, i.e. 
       # mu = x %*% beta
     # If I'm evaluating a stability selection method (which I can tell because
     # the beta is NA), I only have one feature selection set to account for.
     # any(is.na(out$beta)) & 
     if(!is.list(out$selected)){
          return(sc4(out$selected, model$x, proj_x.mu = model$mu))
	} else{
          # Otherwise I am evaluating a glmnet method (lasso or elastic
          # net) and I would like to evaluate SC4 at each lambda,
          # returning the results as a vector (as per the elastic net
          # Simulator vignette). Note that in this case, out$selected
          # is a list of vectors, not a vector.
          list_of_sc4s <- lapply(out$selected, FUN=sc4, X.mat=model$x, 
                                 proj_x.mu=model$mu)
          return(unlist(list_of_sc4s, use.names=FALSE))
     }
     }
)

# Size of selected set

size <- new_metric("size", "Size of Selected Set",
     metric = function(model, out) {
          # # for PCA, simply return the number of principal components
          # if(is.na(out$selected) & !is.na(out$prcomps)){
          #      return(ncol(out$prcomps))
          # }
          # If I'm evaluating a stability selection method (which I can tell
          # because the beta is NA), I only have one feature selection set to
          # account for. Simply find its length.
          if(!is.list(out$selected)){
               return(length(out$selected))
          # Otherwise I am evaluating a glmnet method (lasso or elastic net) and
          # I would like to evaluate the size at each lambda, returning the
          # results as a vector (as per the elastic net Simulator vignette).
          } else{
               return(colSums(as.matrix(out$beta) != 0))
          }
	}
)


# Estimated probabilities

phat <- new_metric("phat", "Estimated Selection Probability",
     metric = function(model, out) {
          if("phat" %in% names(out)){
               return(out$phat)
          } else{
               return(NA)
          }
          
     }
)

# Vector of labels of features

labels <- new_metric("labels", "Feature Labels",
     metric = function(model, out) {
          labels <- rep("Unblocked Noise Features", model$p)
          labels[model$blocked_dgp_vars] <- "Blocked Signal Features"
          labels[model$sig_unblocked_vars] <- "Unblocked Signal Features"
          labels[model$insig_blocked_vars] <- "Blocked Noise Features"

          return(labels)
     }
)

# Whether selected set ontains feature 1

feat1 <- new_metric("feat1", "Feature 1",
     metric = function(model, out) {
          return(as.numeric(1 %in% out$selected))
     }
)

# Whether selected set ontains feature 2

feat2 <- new_metric("feat2", "Feature 2",
     metric = function(model, out) {
          return(as.numeric(2 %in% out$selected))
     }
)

# Whether selected set ontains feature 3

feat3 <- new_metric("feat3", "Feature 3",
     metric = function(model, out) {
          return(as.numeric(3 %in% out$selected))
     }
)

# Whether Equation (16) in notes holds

eq16 <- new_metric("eq16", "Equation 16",
     metric = function(model, out) {
          # Our case imagines X1 was selected first. Since X1 and X2 are
          # exchangeable, for the purpose of this eval we can treat X1
          # as whichever one has a larger correlation with Y.
          cor_1 <- max(out$cor_mat[1, 4], out$cor_mat[2, 4])
          cor_2 <- min(out$cor_mat[1, 4], out$cor_mat[2, 4])
          left_side <- (cor_2 - out$cor_mat[1, 2]*cor_1)/
               (1 - out$cor_mat[1, 2])
          right_side <- (out$cor_mat[3, 4] - out$cor_mat[1, 3]*cor_1)/
               (1 - out$cor_mat[1, 3])
          return(as.numeric(left_side < right_side))
     }
)

# Whether lambda_2^(2) > lambda_1

lambda2_big <- new_metric("lambda2_big", "Lambda2 Big",
     metric = function(model, out) {
          # Our case imagines X1 was selected first. Since X1 and X2 are
          # exchangeable, for the purpose of this eval we can treat X1
          # as whichever one has a larger correlation with Y.
          cor_1 <- max(out$cor_mat[1, 4], out$cor_mat[2, 4])
          cor_2 <- min(out$cor_mat[1, 4], out$cor_mat[2, 4])
          left_side <- (cor_2 - out$cor_mat[1, 2]*cor_1)/
               (1 - out$cor_mat[1, 2])
          return(as.numeric(left_side > cor_1))
     }
)

# Just return correlation matrix

cor_mat_eval <- new_metric("cor_mat_eval", "Correlation Matrix",
     metric = function(model, out) {
          return(out$cor_mat)
     }
)

# Whether event A2 holds

a2 <- new_metric("a2", "A2",
     metric = function(model, out) {
          # Our case imagines X1 was selected first. Since X1 and X2 are
          # exchangeable, for the purpose of this eval we can treat X1
          # as whichever one has a larger correlation with Y.
          cor_1y <- max(out$cor_mat[1, 4], out$cor_mat[2, 4])
          cor_2y <- min(out$cor_mat[1, 4], out$cor_mat[2, 4])
          # We need A1 and S1 to hold, so return true if they don't
          cor_3y <- out$cor_mat[3, 4]
          cor_12 <- out$cor_mat[1, 2]
          if(!(cor_1y > cor_2y & cor_1y > cor_3y & cor_1y > 0 & cor_2y > 0 &
               cor_3y > 0 & cor_12 > 0)){
               return(as.numeric(1))
          }
          return(as.numeric(cor_2y - cor_12*cor_1y >= 0))
     }
)


cssr_mse <- new_metric("cssr_mse", "MSE", metric = function(model, out) {
          out_names <- names(out)

          if("css_res" %in% out_names){
               # This is a method that used the cssr package. We can use
               # getCssPreds to generate predictions.
               return(cssr_mse_metric_func(out, model$sig_blocks +
                    model$k_unblocked))
          } else if("selected_clusts_list" %in% out_names){
               # This is the cluster representative lasso
               return(clus_lasso_metric_func(out, model$sig_blocks +
                    model$k_unblocked))
          } else if("selected_sets" %in% out_names){
               # This is the protolasso
               return(lasso_metric_func(out$selected_sets, out,
                    model$sig_blocks + model$k_unblocked))
          } else if("lasso_selected" %in% out_names){
               # This is the lasso or elastic net
               return(lasso_metric_func(out$lasso_selected, out,
                    model$sig_blocks + model$k_unblocked))
          }

          # Shouldn't be possible to reach this point
          stop("Error: no method for evaluating MSE of this method found")
     }
)

cssr_mse_plant <- new_metric("cssr_mse_plant", "MSE", metric = function(model,
     out) {
          # Get X_train, y_train, X_test, y_test
          data_train <- getPlantSelecData(out$train_inds, model$response_name)
          data_test <- getPlantSelecData(out$test_inds, model$response_name)

          out_names <- names(out)

          if("css_res" %in% out_names){
               # This is a method that used the cssr package. We can use
               # getCssPreds to generate predictions.
               return(cssr_mse_metric_func_plant(out, model,
                    X_train=data_train$X, y_train=data_train$y,
                    X_test=data_test$X, y_test=data_test$y))
          } else if("selected_clusts_list" %in% out_names){
               # This is the cluster representative lasso
               return(clus_lasso_metric_func_plant(out, model,
                    X_train=data_train$X, y_train=data_train$y,
                    X_test=data_test$X, y_test=data_test$y))
          } else if("selected_sets" %in% out_names){
               # This is the protolasso
               return(lasso_metric_func_plant(out$selected_sets, out,
                    model, X_train=data_train$X, y_train=data_train$y,
                    X_test=data_test$X, y_test=data_test$y))
          } else if("lasso_selected" %in% out_names){
               # This is the lasso or elastic net
               return(lasso_metric_func_plant(out$lasso_selected, out,
                    model, X_train=data_train$X, y_train=data_train$y,
                    X_test=data_test$X, y_test=data_test$y))
          }

          # Shouldn't be possible to reach this point
          stop("Error: no method for evaluating MSE of this method found")
     }
)

weight_mse <- new_metric("weight_mse", "Weight MSE", metric = function(model,
     out) {
          out_names <- names(out)

          # This is a setting where all proxies are equally good, so
          # optimal weights are just rep(1/block_size, block_size)

          opt_weights <- rep(1/model$block_size, model$block_size)

          if("css_res" %in% out_names){
               # This is a method that used the cssr package. We can use
               # getCssPreds to generate predictions.
               return(cssr_mse_metric_func(out, model$sig_blocks +
                    model$k_unblocked))
          } else if("selected_clusts_list" %in% out_names){
               # This is the cluster representative lasso; automatically has
               # right weights
               return(0)
          } else if("selected_sets" %in% out_names){
               # This is the protolasso
               return(lasso_metric_func(out$selected_sets, out,
                    model$sig_blocks + model$k_unblocked))
          } else if("lasso_selected" %in% out_names){
               # This is the lasso or elastic net
               return(lasso_metric_func(out$lasso_selected, out,
                    model$sig_blocks + model$k_unblocked))
          }

          # Shouldn't be possible to reach this point
          stop("Error: no method for evaluating MSE of this method found")
     }
)

cssr_mse_metric_func <- function(out, max_model_size){
     # method="ss" # Stability selection
     # "sparse" # sparse cluster stability selection
     # method="weighted_avg" # weighted cluster stability selection
     # method="simple_avg" # simple averaged cluster stability selection

     n_sets <- length(out$selected)
     mses <- rep(as.numeric(NA), max_model_size)
     
     if(n_sets == 0){
          return(mses)
     }

     stopifnot(n_sets <= max_model_size)
     stopifnot(length(out$selected_clusts) >= n_sets)
     stopifnot("selected_clusts" %in% names(out))
     for(i in 1:n_sets){
          # Check if a selected set of size i was defined to exist--if not, skip
          # this model size
          if(!is.null(out$selected[[i]])){
               stopifnot(length(out$selected_clusts[[i]]) == i)
               # Generate test set predictions
               mu_hat_i <- cssr::getCssPreds(out$css_res, testX=out$testX,
                    weighting=out$method, min_num_clusts=i, max_num_clusts=i,
                    trainX=out$testX, trainY=out$testY)
               # Calculate test set MSE
               mses[i] <- mean((mu_hat_i - out$testMu)^2)
          }
     }
     return(mses)
}

cssr_mse_metric_func_plant <- function(out, model, X_train, y_train, X_test,
     y_test){
     # method="ss" # Stability selection
     # "sparse" # sparse cluster stability selection
     # method="weighted_avg" # weighted cluster stability selection
     # method="simple_avg" # simple averaged cluster stability selection

     

     n_sets <- length(out$selected)
     mses <- rep(as.numeric(NA), model$max_model_size)

     if(n_sets == 0){
          return(mses)
     }

     stopifnot(n_sets <= model$max_model_size)
     stopifnot(length(out$selected_clusts) >= n_sets)
     stopifnot("selected_clusts" %in% names(out))
     
     for(i in 1:n_sets){
          # Check if a selected set of size i was defined to exist--if not, skip
          # this model size
          if(!is.null(out$selected[[i]])){
               stopifnot(length(out$selected_clusts[[i]]) == i)
               # Generate test set predictions
               y_hat_i <- cssr::getCssPreds(out$css_res, testX=X_test,
                    weighting=out$method, min_num_clusts=i, max_num_clusts=i,
                    trainX=X_train, trainY=y_train)
               # Calculate test set MSE
               mses[i] <- mean((y_hat_i - y_test)^2)
          }
     }
     return(mses)
}


lasso_metric_func <- function(selected, out, max_model_size){

     n_sets <- max(lengths(selected))
     stopifnot(n_sets <= max_model_size)
     mses <- rep(as.numeric(NA), max_model_size)
     for(i in 1:n_sets){
          inds_i <- which(lengths(selected) == i)
          if(length(inds_i) > 0){
               if(length(inds_i) > 1){
                    sel_set_i <- selected[inds_i][[1]]
               } else{
                    sel_set_i <- selected[[inds_i]]
               }
               mses[i] <- get_mse(out$testX[, sel_set_i], out$testY, out$testMu)
          }
     }
     return(mses)
}

lasso_metric_func_plant <- function(selected, out, model, X_train, y_train,
     X_test, y_test){

     n_sets <- max(lengths(selected))
     stopifnot(n_sets <= model$max_model_size)
     mses <- rep(as.numeric(NA), model$max_model_size)
     for(i in 1:n_sets){
          inds_i <- which(lengths(selected) == i)
          if(length(inds_i) > 0){
               if(length(inds_i) > 1){
                    sel_set_i <- selected[inds_i][[1]]
               } else{
                    sel_set_i <- selected[[inds_i]]
               }
               mses[i] <- get_mse_plant(X_train[, sel_set_i], y_train,
                    X_test[, sel_set_i], y_test)
          }
     }
     return(mses)
}

get_mse <- function(x_train, y_train, mu_train){
     df <- data.frame(y=y_train, x_train)
     df <- df[, colnames(df) != "(Intercept)"]
     lin_model <- stats::lm(y ~. + 0, df)
     preds <- stats::predict.lm(lin_model)
     return(mean((preds - mu_train)^2))
}

get_mse_plant <- function(x_train, y_train, x_test, y_test){
     df_train <- data.frame(y=y_train, x_train)
     df_train <- df_train[, colnames(df_train) != "(Intercept)"]

     df_test <- data.frame(y=y_test, x_test)
     df_test <- df_test[, colnames(df_test) != "(Intercept)"]
     colnames(df_test) <- colnames(df_train)

     lin_model <- stats::lm(y ~., df_train)
     preds <- stats::predict.lm(lin_model, newdata=df_test)
     return(mean((preds - y_test)^2))
}

clus_lasso_metric_func <- function(out, max_model_size){
     n_sets <- length(out$selected_clusts_list)
     stopifnot(n_sets <= max_model_size)
     mses <- rep(as.numeric(NA), max_model_size)
     n <- length(out$testY)
     for(i in 1:n_sets){
          if(!is.null(out$selected_clusts_list[[i]])){
               clusts_list_i <- out$selected_clusts_list[[i]]
               stopifnot(length(clusts_list_i) == i)

               X_train_i <- matrix(as.numeric(NA), nrow=n, ncol=i)

               for(j in 1:i){
                    clust_j_feats <- clusts_list_i[[j]]
                    if(length(clust_j_feats) == 1){
                         X_train_i[, j] <- out$testX[, clust_j_feats]
                    } else if(length(clust_j_feats) > 1){
                         X_train_i[, j] <- rowMeans(out$testX[, clust_j_feats])
                    }
                    
               }
               stopifnot(all(!is.na(X_train_i)))

               mses[i] <- get_mse(X_train_i, out$testY, out$testMu)
          }
     }
     return(mses)
}

clus_lasso_metric_func_plant <- function(out, model, X_train, y_train, X_test,
     y_test){
     n_sets <- length(out$selected_clusts_list)
     stopifnot(n_sets <= model$max_model_size)
     mses <- rep(as.numeric(NA), model$max_model_size)
     n_train <- length(y_train)
     n_test <- length(y_test)
     for(i in 1:n_sets){
          if(!is.null(out$selected_clusts_list[[i]])){
               clusts_list_i <- out$selected_clusts_list[[i]]
               stopifnot(length(clusts_list_i) == i)

               X_train_i <- matrix(as.numeric(NA), nrow=n_train, ncol=i)
               X_test_i <- matrix(as.numeric(NA), nrow=n_test, ncol=i)

               for(j in 1:i){
                    clust_j_feats <- clusts_list_i[[j]]
                    if(length(clust_j_feats) == 1){
                         X_train_i[, j] <- X_train[, clust_j_feats]
                         X_test_i[, j] <- X_test[, clust_j_feats]
                    } else if(length(clust_j_feats) > 1){
                         X_train_i[, j] <- rowMeans(X_train[, clust_j_feats])
                         X_test_i[, j] <- rowMeans(X_test[, clust_j_feats])
                    }   
               }
               stopifnot(all(!is.na(X_train_i)))
               stopifnot(all(!is.na(X_test_i)))

               mses[i] <- get_mse_plant(X_train_i, y_train, X_test_i, y_test)
          }
     }
     return(mses)
}




