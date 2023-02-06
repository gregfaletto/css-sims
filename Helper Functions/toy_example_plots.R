# Going to make plot for slides: get selected set from stability selection,
# first lasso fit

setwd(sim_dir)

set.seed(seed2)

print("n_param_combs:")
print(n_param_combs)

for(p_prime in 1:n_param_combs){
# for(p_prime in 1:1){
    print("p_prime:")
    print(p_prime)

    # Get parameters


    if(n_param_combs > 1){
        output_p_prime <- output(completed_sim)[[p_prime]]
        params_p_prime <- model(completed_sim)[[p_prime]]@params
    }else{
        output_p_prime <- output(completed_sim)
        params_p_prime <- model(completed_sim)@params

    }

    n_model_p_prime <- params_p_prime$n
    p_p_prime <- params_p_prime$p
    k_unblocked_p_prime <- params_p_prime$k_unblocked
    beta_high_p_prime <- params_p_prime$beta_high
    block_size_p_prime <- params_p_prime$block_size
    Sigma_p_prime <- params_p_prime$Sigma

    # Generating plots justifying GSS for paper--skip unneeded plots
    cond1 <- beta_high_p_prime == 1.5 & (p_p_prime %in% c(100, 200, 400))
    cond2 <- beta_high_p_prime == 2 & p_p_prime == 100

    if(!cond1 & !cond2){
        next
    }

    # # For weighted averaging
    # cond3 <- beta_high_p_prime == 2 
    # cond4 <- abs(beta_high_p_prime - sqrt(10)) < 10^(-6)

    # if(!cond3 & !cond4){
    #     next
    # }

    

    if(ncol(Sigma_p_prime) != nrow(Sigma_p_prime)){
        stop("ncol(Sigma_p_prime) != nrow(Sigma_p_prime)")
    }

    if(ncol(Sigma_p_prime) != p_p_prime + 1){
        stop("ncol(Sigma_p_prime) != p_p_prime + 1")
    }

    Sigma_p_prime <- Sigma_p_prime[2:(p_p_prime + 1), 2:(p_p_prime + 1)]





    # Skip this plot if doesen't match p_plot, beta_plot (if specified)

    if(!is.na(p_plot)){
        if(p_p_prime != p_plot){
            next
        }
    }

    if(!is.na(beta_plot)){
        if(beta_high_p_prime != beta_plot){
            next
        }
    }

    # Matrix to store results in

    losses_mat <- matrix(NA, nrow=n_meths_to_eval*p_max, ncol=n_sims)
    losses_mat_rows <- character()
    for(i in 1:n_meths_to_eval){
        loss_name_i <- paste(methods_to_eval[i], "loss")
        losses_mat_rows <- c(losses_mat_rows, paste(loss_name_i, 1:p_max))
    }
    # losses_mat_rows <- c(losses_mat_rows, "True Model loss")
    if(length(losses_mat_rows) != nrow(losses_mat)){
        stop("length(losses_mat_rows) != nrow(losses_mat)")
    }
    rownames(losses_mat) <- losses_mat_rows




    # # Matrix to store cluster counted model sizes in
    # j_loss_cluster_count_mat <- matrix(0, nrow=n_meths_to_eval*p_max,
    #     ncol=n_sims)

    # Matrix to store which model sizes existed for each sim
    j_choices_mat <- matrix(FALSE, nrow=n_meths_to_eval*p_max, ncol=n_sims)

    # Matrix to store number of false selections in

    false_selecs_mat <- matrix(0, nrow=n_meths_to_eval*p_max, ncol=n_sims)

    false_selecs_mat_rows <- character()
    for(i in 1:n_meths_to_eval){
        loss_name_i <- paste("Num", methods_to_eval[i], "False Selections")
        false_selecs_mat_rows <- c(false_selecs_mat_rows, paste(loss_name_i, 1:p_max))
    }
    if(length(false_selecs_mat_rows) != nrow(false_selecs_mat)){
        stop("length(false_selecs_mat_rows) != nrow(false_selecs_mat)")
    }
    rownames(false_selecs_mat) <- false_selecs_mat_rows

    # rownames(false_selecs_mat) <- c(paste("Num Lasso False Selections", 1:p_max),
    #     paste("Num SS False Selections", 1:p_max))

    # Matrix to store number of proxy selections in

    proxy_selecs_mat <- matrix(0, nrow=n_meths_to_eval*p_max, ncol=n_sims)

    proxy_selecs_mat_rows <- character()
    for(i in 1:n_meths_to_eval){
        loss_name_i <- paste("Num", methods_to_eval[i], "Proxy Selections")
        proxy_selecs_mat_rows <- c(proxy_selecs_mat_rows, paste(loss_name_i,
            1:p_max))
    }
    if(length(proxy_selecs_mat_rows) != nrow(proxy_selecs_mat)){
        stop("length(proxy_selecs_mat_rows) != nrow(proxy_selecs_mat)")
    }
    rownames(proxy_selecs_mat) <- proxy_selecs_mat_rows

    # rownames(proxy_selecs_mat) <- c(paste("Num Lasso Proxy Selections", 1:p_max),
    #     paste("Num SS Proxy Selections", 1:p_max))


    # List of lists of Matrices to store selected sets in, as binary vectors 
    # (in order to calculate stability metrics)

    # List where selected clusters are only represented by the most selected
    # feature

    sel_mats <- list()
    for(k in 1:n_meths_to_eval){
        list_k <- list()
        for(j in 1:p_max){
            list_k[[j]] <- matrix(0, n_sims, p_p_prime)
        }
        sel_mats[[k]] <- list_k
    }

    # List where selected clusters are fully represented (i.e., a 1 in the
    # column of every feature in the cluster)

    sel_mats_clusts <- list()
    for(k in 1:n_meths_to_eval){
        list_k <- list()

        if(methods_to_eval[k] %in% c("SS_GSS_random_avg",
            "SS_GSS_random_avg_custom", "SS_GSS_random_avg_unwt",
            "SS_GSS_random_avg_unwt_custom", "BRVZ_avg_unwt")){
            max_j <- p_max + block_size - 1
        } else{
            max_j <- p_max
        }

        for(j in 1:max_j){
            list_k[[j]] <- matrix(0, n_sims, p_p_prime)
        }

        sel_mats_clusts[[k]] <- list_k
    }

    # List like sel_mats_clusts, but j (the model size) is the number of fittted
    # coefficients, not the number of selected features 

    sel_mats_clusts_coefs <- list()
    for(k in 1:n_meths_to_eval){
        list_k <- list()
        for(j in 1:p_max){
            list_k[[j]] <- matrix(0, n_sims, p_p_prime)
        }
        sel_mats_clusts_coefs[[k]] <- list_k
    }

    method_indices <- integer(n_meths_to_eval)

    for(k in 1:n_meths_to_eval){
        method_indices[k] <- identifyMethodIndex(methods_to_eval[k],
            output_p_prime)
    }

    for(i in 1:n_sims){

        # ss_method_indices <- integer(n_stab_meths_to_eval)

        # for(j in 1:n_stab_meths_to_eval){
        #     ss_method_indices[j] <- identifyMethodIndex(stab_meths_to_eval[j])
        # }

        # # How many taus are above 0.5 for each stability selected method
        # p_tau <- integer(n_stab_meths_to_eval)

        # taus <- list()

        # for(k in 1:n_stab_meths_to_eval){
        #     results <- output_p_prime[[ss_method_indices[k]]]@out

        #     taus[[k]] <- results[[i]]$phat

        #     # How many taus are above 0.5?
        #     p_tau[k] <- sum(taus[[k]] > 0.5)
        # }

        


        # Gather stability selected sets

        # Possible values of j
        j_choices <- 1:p_max

        # Index corresponding to this i
        id <- paste("r1.", i, sep="")
        # Selected sets for each method
        selected_sets <- list()
        # If using averaging method, need lists for features to average,
        # clusters to average, and weights
        if(any(c("SS_GSS_random_avg", "SS_GSS_random_avg_custom",
            "SS_GSS_random_avg_unwt", "SS_GSS_random_avg_unwt_custom",
            "BRVZ_avg_unwt") %in% methods_to_eval)){
            to_avg_list_list <- list()
            selected_clusts_list_list <- list()
            weights_list_list <- list()
        }
        # Extract selected sets from each method
        for(k in 1:n_meths_to_eval) {
            # print("methods_to_eval[k]:")
            # print(methods_to_eval[k])
            if(methods_to_eval[k] %in% c("SS_GSS_random_avg",
                "SS_GSS_random_avg_custom", "SS_GSS_random_avg_unwt",
                "SS_GSS_random_avg_unwt_custom", "BRVZ_avg_unwt")){
                selected_results <- extractSelectedSets(
                    output_p_prime[[method_indices[k]]]@out[[id]],
                    methods_to_eval[k], j_choices)
                if(all(is.null(selected_results[[1]]))){
                    stop("all(is.null(selected_results[[1]]))")
                }
                selected_sets[[k]] <- selected_results$selected_sets
                j_choices_k <- selected_results$j_choices

                if(max(j_choices_k) != length(selected_sets[[k]])){
                    print("j_choices_k:")
                    print(j_choices_k)
                    print("max(j_choices_k):")
                    print(max(j_choices_k))
                    print("selected_sets[[k]]:")
                    print(selected_sets[[k]])
                    print("length(selected_sets[[k]]):")
                    print(length(selected_sets[[k]]))
                    stop("max(j_choices_k) != length(selected_sets[[k]])")
                }

                if(sum(lengths(selected_sets[[k]]) > 0) != length(j_choices_k)){
                    stop("sum(lengths(selected_sets[[k]]) > 0) != length(j_choices_k)")
                }

                if(any(lengths(selected_sets[[k]])[j_choices_k] < 1)){
                    stop("any(lengths(selected_sets[[k]])[j_choices_k] < 1)")
                }

                to_avg_list_list[[k]] <- selected_results$to_avg_list
                selected_clusts_list_list[[k]] <- selected_results$selected_clusts_list
                weights_list_list[[k]] <- selected_results$weights_list

                rm(selected_results)

                # Outputs:
                #
                # to_avg_list_list: A list of length n_meths_to_eval. Each
                # element is a list to_avg_list whose elements are logical
                # vectors of length j. For method k, if feature l in 1:j is in a 
                # cluster, the lth entry of the jth vector of the kth
                # to_avg_list in to_avg_list_list will be TRUE. 
                #
                # selected_clusts_list_list: A list of length n_meths_to_eval. Each
                # element is a list selected_clusts_list whose elements are lists of
                # length j. For method k, if feature l in 
                # 1:j is in a cluster, the lth entry of the jth vector of the kth
                # selected_clusts_list in selected_clusts_list_list
                # will contain the features from the
                # cluster of which feature l is a member.
                #
                # weights_list_list: A list of length n_meths_to_eval. Each
                # element is a list weights_list whose elements are lists of
                # length j. For method k, if feature l in 
                # 1:j is in a cluster, the lth entry of the jth vector of the
                # kth weights_list in weights_list_list will
                # be the weights to use (in the same order as the jth entry of
                # the kth selected_clusts_list in selected_clusts_list_list).

            } else{
                # print("(in else)")
                selected_results <- extractSelectedSets(
                    output_p_prime[[method_indices[k]]]@out[[id]],
                    methods_to_eval[k], j_choices)
                if(all(is.null(selected_results[[1]]))){
                    stop("all(is.null(selected_results[[1]]))")
                }
                selected_sets[[k]] <- selected_results$selected_sets
                j_choices_k <- selected_results$j_choices

                rm(selected_results)

                if(max(j_choices_k) != length(selected_sets[[k]])){
                    stop("max(j_choices_k) != length(selected_sets[[k]])")
                }

                if(sum(lengths(selected_sets[[k]]) > 0) != length(j_choices_k)){
                    print("selected_sets[[k]]:")
                    print(selected_sets[[k]])
                    print("lengths(selected_sets[[k]]):")
                    print(lengths(selected_sets[[k]]))
                    print("sum(lengths(selected_sets[[k]]) > 0:")
                    print(sum(lengths(selected_sets[[k]]) > 0))
                    print("j_choices_k:")
                    print(j_choices_k)
                    print("length(j_choices_k)):")
                    print(length(j_choices_k))
                    stop("sum(lengths(selected_sets[[k]]) > 0) != length(j_choices_k)")
                }

                if(any(lengths(selected_sets[[k]])[j_choices_k] < 1)){
                    stop("any(lengths(selected_sets[[k]])[j_choices_k] < 1)")
                }
            }


            # if(methods_to_eval[k]=="SS_GSS_random"){
            #     print(paste("selected_sets for", nameMap(methods_to_eval[k]), ":"))
            #     print(selected_sets[[k]])

            #     print(paste("j_choices for", nameMap(methods_to_eval[k]), ":"))
            #     print(j_choices_k)
            # }


            # print(paste("selected_sets[[", j, "]]:", sep=""))
            # print(selected_sets[[j]])
            # print("selected_sets (in loop):")
            # print(selected_sets)

            # Keep track of these j_choices
            if(any(j_choices_mat[j_choices_k + (k-1)*p_max, i])){
                stop("any(j_choices_mat[j_choices_k + (k-1)*p_max, i])")
            }
            j_choices_mat[j_choices_k + (k-1)*p_max, i] <- TRUE
            
            if(length(j_choices_k) == 0){
                print(paste("length(j_choices_k) == 0 on iteration", i,
                    "for method", methods_to_eval[k]))
                next
            }

            # if(methods_to_eval[k]=="SS_GSS_random"){
            #     print(paste("j_choices_k for:", nameMap(methods_to_eval[k])))
            #     print(j_choices_k)
            # }

            # Save selected sets to later calculate stability metric
            for(j in j_choices_k){
                # if(j <= length(selected_sets[[k]])){
                selected_feat_inds_k_j <- selected_sets[[k]][[j]]

                # if(methods_to_eval[k]=="SS_GSS_random"){
                #     print(paste("j for", nameMap(methods_to_eval[k]), ":"))
                #     print(j)

                #     print(paste("selected_feat_inds_k_j for",
                #         nameMap(methods_to_eval[k]), ":"))
                #     print(selected_feat_inds_k_j)
                # }
            
                if(length(selected_feat_inds_k_j) != j){
                    stop("length(selected_feat_inds_k_j) != j")
                }
                if(any(sel_mats[[k]][[j]][i, ] != 0)){
                    stop("any(sel_mats[[k]][[j]][i, ] != 0)")
                }
                sel_mats[[k]][[j]][i, selected_feat_inds_k_j] <- 1
                # Now deal with cluster representation--for averaging methods,
                # if one of the features is a cluster member, include all
                # cluster members as selected features

                j_clust <- j

                if(methods_to_eval[k] %in% c("SS_GSS_random_avg",
                    "SS_GSS_random_avg_custom", "SS_GSS_random_avg_unwt",
                    "SS_GSS_random_avg_unwt_custom", "BRVZ_avg_unwt")){
                    # Are any of the selected features cluster members?

                    # If none of the selected features are cluster members,
                    # the cluster representation is the same as the non-
                    # cluster representation

                    if(any(to_avg_list_list[[k]][[j]])){
                        # Identify the cluster members
                        clust_feats <- which(to_avg_list_list[[k]][[j]])
                        for(l in clust_feats){
                            selected_feat_inds_k_j <- c(selected_feat_inds_k_j,
                                selected_clusts_list_list[[k]][[j]][[l]])
                            selected_feat_inds_k_j <- unique(selected_feat_inds_k_j)
                        }
                        # This should make j_clust larger
                        if(length(selected_feat_inds_k_j) <= j_clust){
                            stop("length(selected_feat_inds_k_j) <= j_clust")
                        }
                        j_clust <- length(selected_feat_inds_k_j)
                    } 
                } 
                # Selection matrix where j is number of selected features
                if(any(sel_mats_clusts[[k]][[j_clust]][i, ] != 0)){
                    stop("any(sel_mats_clusts[[k]][[j_clust]][i, ] != 0)")
                }
                sel_mats_clusts[[k]][[j_clust]][i, selected_feat_inds_k_j] <- 1

                # Selection matrix where j is number of fitted coefficients
                if(any(sel_mats_clusts_coefs[[k]][[j]][i, ] != 0)){
                    stop("any(sel_mats_clusts_coefs[[k]][[j]][i, ] != 0)")
                }
                sel_mats_clusts_coefs[[k]][[j]][i, selected_feat_inds_k_j] <- 1
                # }
            }

            
        }

        # print("selected_sets (after loop):")
        # print(selected_sets)

        # # Lists of lasso selected features and stability selected features
        # lasso_results <- extractSelectedSets(
        #         output_p_prime[[method_indices[1]]]@out,
        #         methods_to_eval[1], j_choices)

        # lasso_selected_sets <- lasso_results[[1]]

        # j_choices <- lasso_results[[2]]

        # stability_selected_results <- extractSelectedSets(
        #     output_p_prime[[method_indices[2]]]@out,
        #     methods_to_eval[2], j_choices, taus[[1]])
        
        # stability_selected_sets <- stability_selected_results[[1]]
        # j_choices <- stability_selected_results[[2]]


        
      
        # Generate new training and test data

        if(sim=="random_ranking"){
            mu_x_sd_train <- gen_mu_x_sd4_ranking2(n_test, p_p_prime,
                k_unblocked_p_prime,
                beta_low, beta_high_p_prime, nblocks, sig_blocks,
                block_size_p_prime, rho, var, snr, sigma_eps_sq)

            # mu_x_sd_test <- gen_mu_x_sd4_ranking2(n_test, p_p_prime,
            #     k_unblocked_p_prime,
            #     beta_low, beta_high_p_prime, nblocks, sig_blocks,
            #     block_size_p_prime, rho_high, rho_low, var, snr, sigma_eps_sq)
        } else if(sim=="random"){
            mu_x_sd_train <- gen_mu_x_sd4(n_test, p_p_prime, k_unblocked_p_prime,
                beta_low, beta_high_p_prime, nblocks, sig_blocks,
                block_size_p_prime, rho, var, snr, sigma_eps_sq)

            # mu_x_sd_test <- gen_mu_x_sd4(n_test, p_p_prime, k_unblocked_p_prime,
            #     beta_low, beta_high_p_prime, nblocks, sig_blocks,
            #     block_size_p_prime, rho, var, snr, sigma_eps_sq)
        } else if(sim=="random_weighted"){
            mu_x_sd_train <- gen_mu_x_sd4_weighted(n_test, p_p_prime,
                k_unblocked_p_prime, beta_low, beta_high_p_prime, nblocks,
                sig_blocks, block_size_p_prime, n_strong_block_vars, rho_high,
                rho_low, var, snr, sigma_eps_sq)

            # mu_x_sd_test <- gen_mu_x_sd4_weighted(n_test, p_p_prime, k_unblocked_p_prime,
            #     beta_low, beta_high_p_prime, nblocks, sig_blocks,
            #     block_size_p_prime, n_strong_block_vars, rho_high, rho_low, var,
            #     snr, sigma_eps_sq)
            
        } else{
            stop("sim type not random_ranking, random, or random_weighted")
        }
     
        mu_train <- mu_x_sd_train$mu
        x_train <- mu_x_sd_train$x
        sd_train <- mu_x_sd_train$sd
        beta <- mu_x_sd_train$beta
        z_train <- mu_x_sd_train$z

        # mu_test <- mu_x_sd_test[[1]]
        # x_test <- mu_x_sd_test[[2]]
        # sd_test <- mu_x_sd_test[[3]]
        # z_test <- mu_x_sd_test[[6]]

      

        y_train <- mu_train + sd_train * rnorm(n_test)
        # y_test <- mu_test + sd_test * rnorm(n_test)


        # Fit OLS fits for lasso selected set, stability selection selected set,
        # and true model

        # # Vectors of lasso selected features and stability selected features

        # lasso_selected_loss <- numeric(min(p_max, min(p_tau)))
        # # stability_selected_loss <- numeric(min(p_max, min(p_tau)))
        # stability_selected_loss <- matrix(0, min(p_max, min(p_tau)),
        #     n_stab_meths_to_eval)

        # for(j in j_choices){
        #     lasso_selected_loss[j] <- get_loss(x_train, y_train, x_test, y_test,
        #         lasso_selected_sets[[j]])
        #     for(k in 1:n_stab_meths_to_eval){
        #         stability_selected_loss[j, k] <- get_loss(x_train, y_train, x_test,
        #             y_test, stability_selected_sets[[j]][[k]])
        #     }
            
        # }


        for(k in 1:n_meths_to_eval){
            j_choices_k <- which(j_choices_mat[1:p_max + (k-1)*p_max, i])
            if(methods_to_eval[k] %in% c("SS_GSS_random_avg", # no longer used
                "SS_GSS_random_avg_custom", # Weighted averaged cluster
                # stability selection
                "SS_GSS_random_avg_unwt", # no longer used
                "SS_GSS_random_avg_unwt_custom", # Simple averaged cluster
                # stability selection
                "BRVZ_avg_unwt" # Cluster representative lasso
                )){
                for(j in j_choices_k){

                    if(!is.na(losses_mat[(k - 1)*p_max + j, i])){
                        stop("!is.na(losses_mat[(k - 1)*p_max + j, i])")
                    }
                    # losses[j, k] <- get_loss(x_train, y_train, x_test,
                    #     y_test, selected_sets[[k]][[j]])
                    # print("Method:")
                    # print(nameMap(methods_to_eval[k]))
                    # print("selected_sets[[k]][[j]]:")
                    # print(selected_sets[[k]][[j]])
                    # print("to_avg_list_list[[k]][[j]]:")
                    # print(to_avg_list_list[[k]][[j]])
                    losses_mat[(k - 1)*p_max + j, i]  <- get_loss_mu(x_train,
                        y_train, mu_train, selected_sets[[k]][[j]],
                        average=TRUE, to_avg=to_avg_list_list[[k]][[j]],
                        avg_feats=selected_clusts_list_list[[k]][[j]],
                        weights=weights_list_list[[k]][[j]])

                    # # Get actual model size (i.e., number of features selected
                    # # from original feature space, not counting selected
                    # # clusters as only one feature)

                    # model_size <- j

                    # if(any(to_avg_list_list[[k]][[j]])){
                    #     for(l in 1:length(to_avg_list_list[[k]][[j]])){
                    #         if(to_avg_list_list[[k]][[j]][l]){
                    #             model_size <- model_size +
                    #                 length(weights_list_list[[k]][[j]][[l]]) - 1
                    #         }
                    #     }
                    # }

                    # j_loss_cluster_count_mat[(k - 1)*p_max + j, i] <- model_size
                }
            } else{
                for(j in j_choices_k){
                    # losses[j, k] <- get_loss(x_train, y_train, x_test,
                    #     y_test, selected_sets[[k]][[j]])

                    if(!is.na(losses_mat[(k - 1)*p_max + j, i])){
                        stop("!is.na(losses_mat[(k - 1)*p_max + j, i])")
                    }


                    losses_mat[(k - 1)*p_max + j, i] <- get_loss_mu(x_train,
                        y_train, mu_train, selected_sets[[k]][[j]])

                    # j_loss_cluster_count_mat[(k - 1)*p_max + j, i] <- j
                }
            }
            
        }
        # # lasso_loss <- get_loss(x_train, y_train, x_test, y_test, lasso_selected)
        # # print(paste("Lasso loss:", lasso_loss))
        
        # stabs_loss <- get_loss(x_train, y_train, x_test, y_test, stability_selected)
        # print(paste("Stability selection loss (tau = 0.51):", stabs_loss))


        # stabs_top_3_loss <- get_loss(x_train, y_train, x_test, y_test,
        #     stability_selected_top_3)
        # print(paste("Stability selection loss (top 3):", stabs_top_3_loss))
        # print("Stability selected model of size 3:")
        # print(stability_selected_top_3)

        # lasso_loss_3 <- get_loss(x_train, y_train, x_test, y_test, lasso_sets_3)
        # print(paste("Lasso loss (first 3 features):", lasso_loss_3))
        # print("Lasso model of size 3:")
        # print(lasso_sets_3)

        # stabs_top_5_loss <- get_loss(x_train, y_train, x_test, y_test,
        #     stability_selected_top_5)
        # print(paste("Stability selection loss (top 5):", stabs_top_5_loss))
        # print("Stability selected model of size 5:")
        # print(stability_selected_top_5)

        # lasso_loss_5 <- get_loss(x_train, y_train, x_test, y_test, lasso_sets_5)
        # print(paste("Lasso loss (first 5 features):", lasso_loss_5))
        # print("Lasso model of size 5:")
        # print(lasso_sets_5)

        # beta_true <- beta[2:length(beta)]
        # if(length(beta_true) != ncol(x_train)){
        #     stop("length(beta_true) != ncol(x_train)")
        # }
        # if(nblocks==1 & sig_blocks==1){
        #     features_true <- c(1, which(beta_true != 0))
        # } else{
        #     stop("!(nblocks==1 & sig_blocks==1)")
        # }

        # x_train_true <- cbind(z_train, x_train)
        # colnames(x_train_true)[1] <- "z"

        # x_test_true <- cbind(z_test, x_test)
        # colnames(x_test_true)[1] <- "z"

        # true_loss <- get_loss(x_train_true, y_train, x_test_true, y_test,
        #     features_true)

        # if(!is.na(losses_mat[nrow(losses_mat), i])){
        #     stop("!is.na(losses_mat[nrow(losses_mat), i])")
        # }

        # losses_mat[nrow(losses_mat), i] <- true_loss
        # print(paste("True model loss:", true_loss))

        # true_mean_loss <- Metrics::mse(actual=y_test, predicted=mu_test)
        # print(paste("True mean loss:", round(true_mean_loss, 2)))   

        # # Store results in matrix
        # losses_mat[j_choices, i] <- lasso_selected_loss[j_choices]
        # for(k in 1:n_stab_meths_to_eval){
        #     losses_mat[p_max + j_choices, i] <-
        #         stability_selected_loss[j_choices, k]
        # }
        # losses_mat[(n_stab_meths_to_eval + 1)*p_max + 1, i] <- true_loss

        # Now get number of false and proxy selections

        # for(k in 1:n_meths_to_eval){
        #     for(j in j_choices){
        #         false_selecs_mat[j + (k-1)*p_max, i] <-
        #             num_false(selected_sets[[k]][[j]])
        #     }
        # }

        # for(k in 1:n_meths_to_eval){
        #     for(j in j_choices){
        #         proxy_selecs_mat[j + (k-1)*p_max, i] <-
        #             num_proxies(selected_sets[[k]][[j]])
        #     }
        # }

        # If this is the first simulation, create plot of results from just this
        # iteration as a representative example
        if(i == 1){
        	# Visualization of stability selected features
        	if(method=="MB"){
    	        phat_plot <- createPhatPlot(sim_model=model(completed_sim),
                    sim_evals=evals(completed_sim),
    	            method="lassoMB_phat_random", ylab="Proportion of Subsamples"
    	            # , tau1=stability_selected_top_3_tau,
    	            # tau2=stability_selected_top_5_tau
                    )
    	    } else{
                phat_plot <- createPhatPlot(sim_model=model(completed_sim),
                    sim_evals=evals(completed_sim),
                    method="SS_SS_random_custom", ylab="Proportion of Subsamples"
                    # , tau1=stability_selected_top_3_tau,
                    # tau2=stability_selected_top_5_tau
                    )
    	    }

            saveFigure(subdir="Slides", plot=phat_plot, size="slide",
                filename="ss_props.pdf")

    		# Create plot like for final results
            # plot_first_iter <- createLossesPlot(lasso_selected_loss,
            #     stability_selected_loss, true_loss, ss_names=stab_meths_to_eval)

            # print(plot_first_iter)

            # plot_first_iter2 <- createLossesPlot2(as.matrix(losses_mat[, 1]),
            #     methods_to_eval, as.matrix(j_choices_mat[, 1]),
            #     n_model_p_prime, p_p_prime, k_unblocked_p_prime,
            #     beta_high_p_prime, J=p_max)

            # print(plot_first_iter2)
    		

    		# # Also create plot of number of false selections at each model
    		# # size for each model
    		# lasso_num_false <- integer(min(p_max, min(p_tau)))
    		# ss_num_false <- integer(min(p_max, min(p_tau)))
    		# for(j in j_choices){
    	 #        lasso_num_false[j] <- num_false(lasso_selected_sets[[j]])
    	 #        ss_num_false[j] <- num_false(stability_selected_sets[[j]])
    	 #    }


            # false_plot_first_iter2 <- createFalseSelectionsPlot2(
            #     as.matrix(false_selecs_mat[, 1]), methods_to_eval,
            #     as.matrix(j_choices_mat[, 1]), n_model_p_prime, p_p_prime,
            #     k_unblocked_p_prime, beta_high_p_prime, J=p_max)

            # print(false_plot_first_iter2)

            # proxy_plot_first_iter2 <- createProxySelectionsPlot2(
            #     as.matrix(proxy_selecs_mat[, 1]), methods_to_eval,
            #     as.matrix(j_choices_mat[, 1]), n_model_p_prime, p_p_prime,
            #     k_unblocked_p_prime, beta_high_p_prime, J=p_max)

            # print(proxy_plot_first_iter2)
        }


    }
    
    plot_overall2 <- createLossesPlot2(losses_mat, methods_to_eval,
        j_choices_mat, n_model=n_model_p_prime, p=p_p_prime,
        k=k_unblocked_p_prime, beta_high=beta_high_p_prime, legend=mse_legend,
        plot_errors=TRUE, subtitle=FALSE)

    print(plot_overall2)

    saveFigure("MSE vs Num Fitted Coefs", plot_overall2, size="small")

    saveFigure(subdir="Slides", plot=plot_overall2, size="slide",
        filename="gss_loss.pdf")

    # plot_overall2_model_size <- createLossesPlot2(losses_mat, methods_to_eval,
    #     j_choices_mat, n_model_p_prime, p_p_prime, k_unblocked_p_prime,
    #     beta_high_p_prime, legend=mse_legend, xaxis="No. Selected Features",
    #     cluster_count_mat=j_loss_cluster_count_mat)

    # print(plot_overall2_model_size)

    # plot_overall_false2 <- createFalseSelectionsPlot2(false_selecs_mat,
    #     methods_to_eval, j_choices_mat, n_model_p_prime, p_p_prime,
    #     k_unblocked_p_prime, beta_high_p_prime, J=p_max, legend=mse_legend)

    # print(plot_overall_false2)

    # plot_overall_proxy2 <- createProxySelectionsPlot2(proxy_selecs_mat,
    #     methods_to_eval, j_choices_mat, n_model_p_prime, p_p_prime,
    #     k_unblocked_p_prime, beta_high_p_prime, J=p_max, legend=mse_legend)

    # print(plot_overall_proxy2)

    # # grid.arrange(plot_overall_false + theme(legend.position="none"),
    # #     plot_overall_proxy, ncol=2)

    # grid.arrange(plot_overall_false2 + theme(legend.position="none"),
    #     plot_overall_proxy2 + theme(legend.position="none"), ncol=2)


    # # Finally, calculate stability metrics and create plot of those.
    # stab_mets <- matrix(0, p_max, n_meths_to_eval)
    # for(j in 1:p_max){
    #     for(k in 1:n_meths_to_eval){
    #         stab_mets[j, k] <- calcHamStab(sel_mats[[k]][[j]], n_sims)
    #     }
    # }

    # colnames(stab_mets) <- methods_to_eval

    # # print("stab_mets:")
    # # print(stab_mets)

    # plot_stability <- createHamStabPlot(stab_mets, n_model_p_prime, p_p_prime, 
    #     k_unblocked_p_prime, beta_high_p_prime)

    # print(plot_stability)


    ########## Stability metrics vs num. fitted coefficients

    # Finally, calculate stability metrics and create plot of those.
    stab_mets_NSB <- matrix(as.numeric(NA), p_max, n_meths_to_eval)
    stab_mets_NSB_conf_lower <- matrix(as.numeric(NA), p_max,
        n_meths_to_eval)
    stab_mets_NSB_conf_upper <- matrix(as.numeric(NA), p_max,
        n_meths_to_eval)

    if(calc_clust_stab){
        # Additional matrices for modified stability metric where all clustered
        # features count the same

        stab_mets_NSB_clust_multi <- matrix(as.numeric(NA), p_max, n_meths_to_eval)
        stab_mets_NSB_conf_lower_clust_multi <- matrix(as.numeric(NA), p_max,
            n_meths_to_eval)
        stab_mets_NSB_conf_upper_clust_multi <- matrix(as.numeric(NA), p_max,
            n_meths_to_eval)

        # Binary indicator for whether cluster was selected

        stab_mets_NSB_clust_single <- matrix(as.numeric(NA), p_max, n_meths_to_eval)
        stab_mets_NSB_conf_lower_clust_single <- matrix(as.numeric(NA), p_max,
            n_meths_to_eval)
        stab_mets_NSB_conf_upper_clust_single <- matrix(as.numeric(NA), p_max,
            n_meths_to_eval)

        # Method proposed by NSB et al in 2020 paper

        stab_mets_NSB_clust_corr <- matrix(as.numeric(NA), p_max, n_meths_to_eval)
        # stab_mets_NSB_conf_lower_clust_corr <- matrix(as.numeric(NA), p_max,
        #     n_meths_to_eval)
        # stab_mets_NSB_conf_upper_clust_corr <- matrix(as.numeric(NA), p_max,
        #     n_meths_to_eval)

    }


    for(j in 1:p_max){
        for(k in 1:n_meths_to_eval){
            if(methods_to_eval[k] %in% c("SS_GSS_random_avg",
                "SS_GSS_random_avg_custom", "SS_GSS_random_avg_unwt",
                "SS_GSS_random_avg_unwt_custom", "BRVZ_avg_unwt")){
                stat_results <- calcNSBStab(sel_mats_clusts_coefs[[k]][[j]],
                    n_sims, calc_errors=TRUE)
            } else{
                stat_results <- calcNSBStab(sel_mats[[k]][[j]], n_sims,
                    calc_errors=TRUE)
            }

            if(!is.na(stab_mets_NSB[j, k])){
                stop("!is.na(stab_mets_NSB[j, k])")
            }

            stab_mets_NSB[j, k] <- stat_results[1]
            stab_mets_NSB_conf_lower[j, k] <- stat_results[2]
            stab_mets_NSB_conf_upper[j, k] <- stat_results[3]

            if(calc_clust_stab){
                if(methods_to_eval[k] %in% c("SS_GSS_random_avg",
                    "SS_GSS_random_avg_custom", "SS_GSS_random_avg_unwt",
                    "SS_GSS_random_avg_unwt_custom", "BRVZ_avg_unwt")){
                    stat_results_clust_multi <- calcNSBStab(sel_mats_clusts_coefs[[k]][[j]],
                        n_sims, calc_errors=TRUE, cluster_count="multi", 
                        n_clust_feats=block_size)

                    stat_results_clust_single <- calcNSBStab(sel_mats_clusts_coefs[[k]][[j]],
                        n_sims, calc_errors=TRUE, cluster_count="single", 
                        n_clust_feats=block_size)

                    stat_results_clust_corr <- calcNSBStab(sel_mats_clusts_coefs[[k]][[j]],
                        n_sims, calc_errors=TRUE, cluster_count="corr", 
                        n_clust_feats=block_size, C=abs(Sigma_p_prime))
                } else{
                    stat_results_clust_multi <- calcNSBStab(sel_mats[[k]][[j]],
                        n_sims, calc_errors=TRUE, cluster_count="multi", 
                        n_clust_feats=block_size)

                    stat_results_clust_single <- calcNSBStab(sel_mats[[k]][[j]],
                        n_sims, calc_errors=TRUE, cluster_count="single", 
                        n_clust_feats=block_size)

                    stat_results_clust_corr <- calcNSBStab(sel_mats[[k]][[j]],
                        n_sims, calc_errors=TRUE, cluster_count="corr", 
                        n_clust_feats=block_size, C=abs(Sigma_p_prime))

                }

                if(!is.na(stab_mets_NSB_clust_multi[j, k])){
                    stop("!is.na(stab_mets_NSB_clust_multi[j, k])")
                }
                
                stab_mets_NSB_clust_multi[j, k] <- stat_results_clust_multi[1]
                stab_mets_NSB_conf_lower_clust_multi[j, k] <- stat_results_clust_multi[2]
                stab_mets_NSB_conf_upper_clust_multi[j, k] <- stat_results_clust_multi[3]

                if(!is.na(stab_mets_NSB_clust_single[j, k])){
                    stop("!is.na(stab_mets_NSB_clust_single[j, k])")
                }
                
                stab_mets_NSB_clust_single[j, k] <- stat_results_clust_single[1]
                stab_mets_NSB_conf_lower_clust_single[j, k] <- stat_results_clust_single[2]
                stab_mets_NSB_conf_upper_clust_single[j, k] <- stat_results_clust_single[3]

                if(!is.na(stab_mets_NSB_clust_corr[j, k])){
                    stop("!is.na(stab_mets_NSB_clust_corr[j, k])")
                }
                
                stab_mets_NSB_clust_corr[j, k] <- stat_results_clust_corr[1]
                # stab_mets_NSB_conf_lower_clust_corr[j, k] <- stat_results_clust_corr[2]
                # stab_mets_NSB_conf_upper_clust_corr[j, k] <- stat_results_clust_corr[3]
            }
        }
    }

    colnames(stab_mets_NSB) <- methods_to_eval

    # print("stab_mets_NSB:")
    # print(stab_mets_NSB)

    plot_stability_NSB_ints <- createNSBStabPlot(stab_mets_NSB,
        n_model=n_model_p_prime, p=p_p_prime, k=k_unblocked_p_prime,
        beta_high=beta_high_p_prime, plot_errors = TRUE,
        lowers=stab_mets_NSB_conf_lower,
        uppers=stab_mets_NSB_conf_upper, subtitle=FALSE)

    print(plot_stability_NSB_ints)

    saveFigure("Stability vs Num Fitted Coefficients/Original Stabilty Metric/Confidence Intervals",
        plot_stability_NSB_ints, size="small")

    saveFigure(subdir="Slides", plot=plot_stability_NSB_ints, size="slide",
        filename="gss_stab.pdf")

    # Side-by-side plots

    # 2. Save the legend
    #+++++++++++++++++++++++
    legend <- get_legend(plot_overall2 + theme(legend.direction="horizontal"))
    # 3. Remove the legend from the box plot
    #+++++++++++++++++++++++
    plot_overall2_no_legend <- plot_overall2 + theme(legend.position="none")
    plot_stability_NSB_ints_no_legend <- plot_stability_NSB_ints +
        theme(legend.position="none")
    blankPlot <- ggplot() + geom_blank(aes(1,1)) + cowplot::theme_nothing()
    # 4. Arrange ggplot2 graphs with a specific width

    side_plot <- grid.arrange(plot_overall2_no_legend,
        plot_stability_NSB_ints_no_legend, legend
        , ncol=3
        , widths=c(0.4, 0.4, 0.2)
        )

    side_plot <- cowplot::ggdraw(side_plot) +
        theme(plot.background = element_rect(fill="white", color = NA))

    print(side_plot)

    saveFigure("MSE and Stability vs Num Fitted Coefficients/Original Stabilty Metric/Confidence Intervals",
        side_plot, size="large")


    # 2. Save the legend
    #+++++++++++++++++++++++
    legend_phat_plot <- get_legend(phat_plot +
        theme(legend.direction="horizontal"))

    legend_plot_overall2 <- get_legend(plot_overall2 +
        theme(legend.direction="horizontal"))
    # 3. Remove the legend from the box plot
    #+++++++++++++++++++++++
    plot_overall2_no_legend <- plot_overall2 + theme(legend.position="none")
    phat_plot_no_legend <- phat_plot + theme(legend.position="none")

    # 4. Arrange ggplot2 graphs with a specific width

    side_plot_intro <- grid.arrange(phat_plot_no_legend,
        plot_overall2_no_legend, legend_phat_plot, legend_plot_overall2
        , ncol=2, nrow=2, widths=c(0.5, 0.5), heights = c(2.5, 0.2)
        )



    # side_plot <- grid.arrange(plot_overall2_no_legend,
    #     plot_stability_NSB_ints_no_legend, plot_stability_NSB_vs_mse_no_legend,
    #     legend, ncol=3, nrow = 2, layout_matrix = rbind(c(1, 2, 3), c(4, 4, 4)),
    #     widths = c(1.8, 1.8, 1.8), heights = c(2.5, 0.2))

    

    #  # 4. Arrange ggplot2 graphs with a specific width

    

    side_plot_intro <- cowplot::ggdraw(side_plot_intro) +
        theme(plot.background = element_rect(fill="white", color = NA))

    print(side_plot)

    saveFigure("Slides", side_plot_intro, size="mlarge")


    plot_stability_NSB <- createNSBStabPlot(stab_mets_NSB,
        n_model_p_prime, p_p_prime, k_unblocked_p_prime, beta_high_p_prime,
        plot_errors = FALSE)

    print(plot_stability_NSB)

    saveFigure("Stability vs Num Fitted Coefficients/Original Stabilty Metric",
        plot_stability_NSB)



    plot_stability_NSB_vs_mse <- createStabMSEPlot(stab_mets_NSB,
        losses_mat, methods_to_eval, j_choices_mat, n_model=n_model_p_prime,
        p=p_p_prime, k=k_unblocked_p_prime, beta_high=beta_high_p_prime,
        plot_errors = FALSE
        # , lowers=stab_mets_NSB_conf_lower
        # , uppers=stab_mets_NSB_conf_upper
        , legend=mse_legend, subtitle=FALSE)

    print(plot_stability_NSB_vs_mse)

    saveFigure("Stability vs MSE/Original Stabilty Metric",
        plot_stability_NSB_vs_mse, size="medium")

    saveFigure(subdir="Slides", plot=plot_stability_NSB_vs_mse, size="slide",
        filename="gss_mse_vs_stab.pdf")


    # 3 figure Side-by-side plot

    # 3. Remove the legend from the box plot
    #+++++++++++++++++++++++
    plot_stability_NSB_vs_mse_no_legend <- plot_stability_NSB_vs_mse +
        theme(legend.position="none")

    # 4. Arrange ggplot2 graphs with a specific width

    side_plot <- grid.arrange(plot_overall2_no_legend,
        plot_stability_NSB_ints_no_legend, plot_stability_NSB_vs_mse_no_legend,
        legend, ncol=3, nrow = 2, layout_matrix = rbind(c(1, 2, 3), c(4, 4, 4)),
        widths = c(1.8, 1.8, 1.8), heights = c(2.5, 0.2))

    side_plot <- cowplot::ggdraw(side_plot) +
        theme(plot.background = element_rect(fill="white", color = NA))

    print(side_plot)

    saveFigure("MSE Stability vs Num Coefs, MSE vs Stability/Original Stabilty Metric/Confidence Intervals",
        side_plot, size="large")





    if(calc_clust_stab){
        # colnames(stab_mets_NSB_clust_multi) <- methods_to_eval

        # # print("stab_mets_NSB:")s
        # # print(stab_mets_NSB)

        # plot_stability_NSB_ints_clust_multi <- createNSBStabPlot(stab_mets_NSB_clust_multi,
        #     n_model_p_prime, p_p_prime, k_unblocked_p_prime, beta_high_p_prime,
        #     plot_errors = TRUE, lowers=stab_mets_NSB_conf_lower_clust_multi,
        #     uppers=stab_mets_NSB_conf_upper_clust_multi, cluster_count="multi")

        # print(plot_stability_NSB_ints_clust_multi)

        # saveFigure("Stability vs Num Fitted Coefficients/Modified Stability (Cluster Counting, multi)/Confidence Intervals",
        #     plot_stability_NSB_ints_clust_multi)

        # plot_stability_NSB_clust_multi <- createNSBStabPlot(stab_mets_NSB_clust_multi,
        #     n_model_p_prime, p_p_prime, k_unblocked_p_prime, beta_high_p_prime,
        #     plot_errors = FALSE, cluster_count="multi")

        # print(plot_stability_NSB_clust_multi)

        # saveFigure("Stability vs Num Fitted Coefficients/Modified Stability (Cluster Counting, multi)",
        #     plot_stability_NSB_clust_multi)


        # plot_stability_NSB_clust_multi_vs_mse <- createStabMSEPlot(stab_mets_NSB_clust_multi,
        #     losses_mat, methods_to_eval, j_choices_mat, n_model_p_prime, p_p_prime,
        #     k_unblocked_p_prime, beta_high_p_prime, plot_errors = FALSE,
        #     legend=mse_legend, cluster_count="multi")

        # print(plot_stability_NSB_clust_multi_vs_mse)

        # saveFigure("Stability vs MSE/Modified Stability (Cluster Counting, multi)",
        #     plot_stability_NSB_clust_multi_vs_mse)







        # colnames(stab_mets_NSB_clust_single) <- methods_to_eval

        # # print("stab_mets_NSB:")s
        # # print(stab_mets_NSB)

        # plot_stability_NSB_ints_clust_single <- createNSBStabPlot(stab_mets_NSB_clust_single,
        #     n_model_p_prime, p_p_prime, k_unblocked_p_prime, beta_high_p_prime,
        #     plot_errors = TRUE, lowers=stab_mets_NSB_conf_lower_clust_single,
        #     uppers=stab_mets_NSB_conf_upper_clust_single, cluster_count="single")

        # print(plot_stability_NSB_ints_clust_single)

        # saveFigure("Stability vs Num Fitted Coefficients/Modified Stability (Cluster Counting, single)/Confidence Intervals",
        #     plot_stability_NSB_ints_clust_single)

        # plot_stability_NSB_clust_single <- createNSBStabPlot(stab_mets_NSB_clust_single,
        #     n_model_p_prime, p_p_prime, k_unblocked_p_prime, beta_high_p_prime,
        #     plot_errors = FALSE, cluster_count="single")

        # print(plot_stability_NSB_clust_single)

        # saveFigure("Stability vs Num Fitted Coefficients/Modified Stability (Cluster Counting, single)",
        #     plot_stability_NSB_clust_single)

        # plot_stability_NSB_clust_single_vs_mse <- createStabMSEPlot(stab_mets_NSB_clust_single,
        #     losses_mat, methods_to_eval, j_choices_mat, n_model_p_prime, p_p_prime,
        #     k_unblocked_p_prime, beta_high_p_prime, plot_errors = FALSE,
        #     legend=mse_legend, cluster_count="single")

        # print(plot_stability_NSB_clust_single_vs_mse)

        # saveFigure("Stability vs MSE/Modified Stability (Cluster Counting, single)",
        #     plot_stability_NSB_clust_single_vs_mse)







        colnames(stab_mets_NSB_clust_corr) <- methods_to_eval

        # print("stab_mets_NSB:")s
        # print(stab_mets_NSB)

        plot_stability_NSB_clust_corr <- createNSBStabPlot(
            stab_mets_NSB_clust_corr,
            n_model=n_model_p_prime, p=p_p_prime, k=k_unblocked_p_prime, 
            beta_high=beta_high_p_prime,
            plot_errors = FALSE, cluster_count="corr", subtitle=FALSE)

        print(plot_stability_NSB_clust_corr)

        saveFigure("Stability vs Num Fitted Coefficients/2020 Correlated Stability",
            plot_stability_NSB_clust_corr, size="small")


        plot_stability_NSB_clust_corr_vs_mse <- createStabMSEPlot(stab_mets_NSB_clust_corr,
            losses_mat, methods_to_eval, j_choices_mat, n_model=n_model_p_prime,
            p=p_p_prime, k=k_unblocked_p_prime, beta_high=beta_high_p_prime,
            plot_errors = TRUE, legend=mse_legend, cluster_count="corr",
            subtitle=FALSE)

        print(plot_stability_NSB_clust_corr_vs_mse)

        saveFigure("Stability vs MSE/2020 Correlated Stability",
            plot_stability_NSB_clust_corr_vs_mse, size="medium")

    }


    


    # ################### Stability metrics vs num. selected features

    # if(any(c("SS_GSS_random_avg", "SS_GSS_random_avg_unwt",
    #     "BRVZ_avg_unwt") %in% methods_to_eval)){
    #     # Finally, calculate stability metrics and create plot of those.
    #     stab_mets_NSB_full <- matrix(as.numeric(NA), p_max + block_size - 1, n_meths_to_eval)
    #     stab_mets_NSB_conf_lower_full <- matrix(as.numeric(NA), p_max + block_size - 1,
    #         n_meths_to_eval)
    #     stab_mets_NSB_conf_upper_full <- matrix(as.numeric(NA), p_max + block_size - 1,
    #         n_meths_to_eval)

    #     if(calc_clust_stab){
    #         # Additional matrices for modified stability metric where all clustered
    #         # features count the same

    #         stab_mets_NSB_clust_full <- matrix(as.numeric(NA), p_max + block_size - 1, n_meths_to_eval)
    #         stab_mets_NSB_conf_lower_clust_full <- matrix(as.numeric(NA), p_max + block_size - 1,
    #             n_meths_to_eval)
    #         stab_mets_NSB_conf_upper_clust_full <- matrix(as.numeric(NA), p_max + block_size - 1,
    #             n_meths_to_eval)
    #     }


    #     for(j in 1:(p_max + block_size - 1)){
    #         for(k in 1:n_meths_to_eval){
    #             if((j <= p_max) | (methods_to_eval[k] %in% c("SS_GSS_random_avg",
    #             "SS_GSS_random_avg_unwt", "BRVZ_avg_unwt"))){
    #                 if(length(sel_mats_clusts[[k]][[j]]) > 0){
    #                     stat_results <- calcNSBStab(sel_mats_clusts[[k]][[j]],
    #                         n_sims, calc_errors=TRUE)

    #                     if(!is.na(stab_mets_NSB_full[j, k])){
    #                         stop("!is.na(stab_mets_NSB_full[j, k])")
    #                     }
    #                     stab_mets_NSB_full[j, k] <- stat_results[1]
    #                     stab_mets_NSB_conf_lower_full[j, k] <- stat_results[2]
    #                     stab_mets_NSB_conf_upper_full[j, k] <- stat_results[3]

    #                     if(calc_clust_stab){
    #                         stat_results_clust <- calcNSBStab(sel_mats_clusts[[k]][[j]],
    #                             n_sims, calc_errors=TRUE, cluster_count="multi", 
    #                             n_clust_feats=block_size)

    #                         if(!is.na(stab_mets_NSB_clust_full[j, k])){
    #                             stop("!is.na(stab_mets_NSB_clust_full[j, k])")
    #                         }

    #                         stab_mets_NSB_clust_full[j, k] <- stat_results_clust[1]
    #                         stab_mets_NSB_conf_lower_clust_full[j, k] <- stat_results_clust[2]
    #                         stab_mets_NSB_conf_upper_clust_full[j, k] <- stat_results_clust[3]
    #                     }
    #                 }
    #             }
    #         }
    #     }

    #     colnames(stab_mets_NSB_full) <- methods_to_eval

    #     # print("stab_mets_NSB:")
    #     # print(stab_mets_NSB)

    #     plot_stability_NSB_ints_full <- createNSBStabPlot(stab_mets_NSB_full,
    #         n_model_p_prime, p_p_prime, k_unblocked_p_prime, beta_high_p_prime,
    #         plot_errors = TRUE, lowers=stab_mets_NSB_conf_lower_full,
    #         uppers=stab_mets_NSB_conf_upper_full, xaxis="No. Selected Features")

    #     print(plot_stability_NSB_ints_full)

    #     plot_stability_NSB_full <- createNSBStabPlot(stab_mets_NSB_full,
    #         n_model_p_prime, p_p_prime, k_unblocked_p_prime, beta_high_p_prime,
    #         plot_errors = FALSE, xaxis="No. Selected Features")

    #     print(plot_stability_NSB_full)

    #     if(calc_clust_stab){
    #         colnames(stab_mets_NSB_clust_full) <- methods_to_eval

    #         # print("stab_mets_NSB:")s
    #         # print(stab_mets_NSB)

    #         plot_stability_NSB_ints_clust_full <- createNSBStabPlot(stab_mets_NSB_clust_full,
    #             n_model_p_prime, p_p_prime, k_unblocked_p_prime, beta_high_p_prime,
    #             plot_errors = TRUE, lowers=stab_mets_NSB_conf_lower_clust_full,
    #             uppers=stab_mets_NSB_conf_upper_clust_full, 
    #             cluster_count="multi", xaxis="No. Selected Features")

    #         print(plot_stability_NSB_ints_clust_full)

    #         plot_stability_NSB_clust_full <- createNSBStabPlot(stab_mets_NSB_clust_full,
    #             n_model_p_prime, p_p_prime, k_unblocked_p_prime, beta_high_p_prime,
    #             plot_errors = FALSE, 
    #             cluster_count="multi", xaxis="No. Selected Features")

    #         print(plot_stability_NSB_clust_full)

        # }

    # }

}

    

    # get out-of-sample predictions and get out-of-sample error.


    # plotPredictionErrors(sim=gss_random_ranking, method="SS_SS_random")

    # plotPredictionErrors(sim=gss_random_ranking, method="SS_GSS_random")

    # plotPredictionErrors(sim=gss_random_ranking, method="lassoSS_phat_cor")

    # #### Find corresponding prediction errors

    # # Add a column to data.frame for prediction errors
    # colnames_orig <- colnames(data)
    # data <- data.frame(data, numeric(p), numeric(p), numeric(p), numeric(p),
    #     numeric(p))
    # colnames(data) <- c(colnames_orig, "PredictionErrorMu",
    #     "ScaledPredictionErrorMu", "PredictionErrorY", "ScaledPredictionErrorY",
    #     "ScreeningCriterion")

    # colnames_orig_ideal <- colnames(data_ideal)
    # data_ideal <- data.frame(data_ideal, numeric(p), numeric(p), numeric(p), numeric(p),
    #     numeric(p))
    # colnames(data_ideal) <- c(colnames_orig_ideal, "PredictionErrorMu",
    #     "ScaledPredictionErrorMu", "PredictionErrorY", "ScaledPredictionErrorY",
    #     "ScreeningCriterion")

    # # True model
    # mu <- model(gss_random_ranking)@params$mu

    # # Observed responses
    # Y <- draws(gss_random_ranking)@draws$r1.1

    # # For k in 1:p, calculate prediction error for top k variables as ranked by
    # # stability selection relative to true model.
    # for(k in 1:p){
    #     # Identify indices of top k features
    #     features_k <- data$Feature[1:k]
    #     X <- as.matrix(model(gss_random_ranking)@params$x[, features_k])
    #     # Predictions from OLS based on top k variables as ranked by stability
    #     # selection, using mu as response
    #     predictions_mu_k <- lm.fit(x=X, y=mu, tol=1e-12)$fitted.values
    #     # Calculate MSE
    #     error_mu_k <- (1/n)*((predictions_mu_k - mu)%*%(predictions_mu_k - mu))
    #     data[k, "PredictionErrorMu"] <- error_mu_k
    #     data[k, "ScaledPredictionErrorMu"] <- error_mu_k/max(data$PredictionErrorMu)
    #     # Predictions from OLS based on top k variables as ranked by stability
    #     # selection, using y as response
    #     predictions_y_k <- lm.fit(x=X, y=Y, tol=1e-12)$fitted.values
    #     # Calculate MSE
    #     error_y_k <- (1/n)*((predictions_y_k - Y)%*%(predictions_y_k - Y))
    #     data[k, "PredictionErrorY"] <- error_y_k
    #     data[k, "ScaledPredictionErrorY"] <- error_y_k/max(data$PredictionErrorY)

    #     # Screening criterion of top k variables as ranked by stability
    #     # selection
    #     sc4_k <- sum(predictions_mu_k^2)/sum(mu^2)
    #     data[k, "ScreeningCriterion"] <- sc4_k
    # }

    # # plot2a <- ggplot(data, aes(x=Ranking, y=PredictionErrorMu)) + geom_line()
    # # # print(plot2a)

    # # plot3 <- ggplot(data, aes(x=Ranking, y=PredictionErrorY)) + geom_line()
    # # print(plot3)

    # plot4 <- ggplot(data) + geom_col(mapping=aes(x=Ranking, y=phat, fill=Label)) +
    #    geom_line(mapping=aes(x=Ranking, y=ScaledPredictionErrorY))
    # print(plot4)

    # # plot4b <- ggplot(data_ideal, aes(x=Ranking, y=phat, fill=Label)) + geom_col() +
    # #    geom_line(aes(x=Ranking, y=ScaledPredictionErrorY))
    # # print(plot4b)

    # plot5 <- ggplot(data) + geom_col(mapping=aes(x=Ranking, y=phat, fill=Label)) +
    #    geom_line(mapping=aes(x=Ranking, y=ScreeningCriterion))
    # print(plot5)

    # # plot5b <- ggplot(data_ideal, aes(x=Ranking, y=phat, fill=Label)) + geom_col() +
    # #    geom_line(aes(x=Ranking, y=ScreeningCriterion))
    # # print(plot5b)

if(sim == "random_ranking" & "SS_GSS_random_custom" %in% methods_to_eval){
    # Print out unusual results from sparse CSS model of size 1 in ranking
    # simulation

    sel_mat <- sel_mats[[which(methods_to_eval=="SS_GSS_random_custom")]][[1]]
    # Eliminate rows of all zeroes (no selections)

    all_zero_rows <- apply(sel_mat, 1, function(x){all(x == 0)})
    print(paste("Selected sets for", nameMap("SS_GSS_random_custom"),
        "model of size 1:"))
    print(sel_mat[!all_zero_rows, ])
}

setwd(wd)

print("Complete! Total time:")
print(Sys.time() - t0)
