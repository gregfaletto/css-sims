############# Utility Functions

complementaryPairSamplingIndicators <- function(X.dat, y, B, q){
	# Output: returns a p x 2B matrix of selected features. Each column contains
	# a set of selected features, with a 1 in row i of column j if feature i was
	# selected in iteration j, and a 0 otherwise.
	#
	# Inputs:
	# X.dat: design matrix
	# y: responses
	# B: number of resamples
	# q: Number of features to select on each lasso iteration (glmnet chooses a 
	# suitable lambda such that (approximately) q features are selected)

	if(q != round(q) | q <= 0){
		stop("q != round(q) | q <= 0")
	}
	n <- nrow(X.dat)
	p <- ncol(X.dat)
	# Per Samworth & Shah, take a sample (without replacement)
	# of size floor(n/2)
	samp_size <- floor(n/2)
	# p by 2B matrix; in each of 2B columns, the jth row will equal 1
	# if that feature was selected, and zero otherwise
	selected_features <- matrix(numeric(p*2*B), ncol=(2*B))
	# Complete B repetitions of the following:
	for(i in 1:B){
		# Per Samworth & Shah, take complementary pairwise samples
		order <- sample(n)
		sample1 <- order[1:samp_size]
		sample2 <- order[(samp_size+1):(samp_size+samp_size)]
		# Fit the lasso using glmnet, choosing range of lambdas such that a
		# maximum of q features are selected.
		# use suppressWarnings because otherwise it gives this useless
		# error about q working properly
		q_to_check <- c(q, q + 1, q - 1, q + 2, q - 2, q + 3, q - 3)
		q_to_check <- q_to_check[q_to_check > 0]

		selected1_q <- numeric()
		nlambda <- 1000*q
		lambda.min.ratio <- max(0.01, 1/(2*q^2))
		nloops <- 1

		selected2_q <- numeric()
		nlambda <- 1000*q
		# lambda.min.ratio <- max(0.01, 1/(2*q^1.5))
		# nloops <- 1

		while(length(selected1_q) == 0 | length(selected2_q) == 0){
			fit1 <- suppressWarnings(glmnet(x=X.dat[sample1, ], y=y[sample1],
				family="gaussian", alpha=1, nlambda=nlambda, lambda.min.ratio=lambda.min.ratio))
			selected1 <- unique(predict(fit1, type="nonzero"))
			# Lasso selected set of size q
	        selected1_q <- selected1[lapply(selected1, length) %in% q_to_check]


	        fit2 <- suppressWarnings(glmnet(x=X.dat[sample2, ], y=y[sample2],
				family="gaussian", alpha=1, nlambda=1000, lambda.min.ratio=0.25))
			selected2 <- unique(predict(fit2, type="nonzero"))
			# Lasso selected set of size q
	        selected2_q <- selected2[lapply(selected2, length) %in% q_to_check]
	        # if(length(selected1_q) == 0){
	        # 	print(unique(predict(fit1, type="nonzero")))
	        # 	print("q:")
	        # 	print(q)
	        #     stop("length(selected1_q) == 0")
	        # }
	        # If we didn't get a match, try again with a new sample, bigger
	        # nlambda, and smaller lambda.min.ratio.
	        if(nloops < 4){
	        	nlambda <- nlambda*10
	        	lambda.min.ratio <- lambda.min.ratio/3
	        }
	        if(nloops > 2){
	        	print("nloops:")
	        	print(nloops)
	        }
	        nloops <- nloops + 1
	        order <- sample(n)
			sample1 <- order[1:samp_size]
			sample2 <- order[(samp_size+1):(samp_size+samp_size)]
	        if(nloops == 10){
	        	print(X.dat)
	        	stop("nloops too big")
	        }
        }
        for(q_prime in q_to_check){
        	if(any(lapply(selected1, length) == q_prime)){
        		set_to_choose <- which(lapply(selected1, length) == q_prime)
        		# print("set_to_choose:")
        		# print(set_to_choose)
				set_to_choose <- set_to_choose[1]
				selected1 <- selected1[[set_to_choose]]
				break
        	}
        	if(length(selected1) == 0){
        		next
        	}
        }
        if(length(selected1) == 0){
        	print(unique(predict(fit1, type="nonzero")))
        	print("q:")
        	print(q)
    		stop("length(selected2) == 0")
    	}

    	

        for(q_prime in q_to_check){
        	if(any(lapply(selected2, length) == q_prime)){
        		set_to_choose <- which(lapply(selected2, length) == q_prime)
				set_to_choose <- set_to_choose[1]
				selected2 <- selected2[[set_to_choose]]
				break
        	}
        	if(length(selected2) == 0){
        		next
        	}
        }
        if(length(selected2) == 0){
        	print(unique(predict(fit2, type="nonzero")))
        	print("q:")
        	print(q)
    		stop("length(selected2) == 0")
    	}

		# fit1 <- suppressWarnings(glmnet(x=X.dat[sample1, ],
		# 	y=y[sample1], family="gaussian", pmax=q))
		# fit2 <- suppressWarnings(glmnet(x=X.dat[sample2, ],
		# 	y=y[sample2], family="gaussian", pmax=q))
		# ## which coefficients are non-zero?
	 #    selected1 <- predict(fit1, type="nonzero")
	 #    selected2 <- predict(fit2, type="nonzero")
	 #    # Choose the features chosen at the particular lambda where q
	 #    # features were selected (i.e. the smallest or last lambda)
	 #    for(i in length(selected1):1){
	 #    	if(length(selected1[[i]]) != 0){
	 #    		selected1 <- selected1[[i]]
	 #    		# print(selected1)
	 #    		break
	 #    	}
	 #    	if(i==1){
	 #    		print(selected1)
	 #    		stop("All selected sets empty")
	 #    	}
	 #    }
	 #    for(i in length(selected2):1){
	 #    	if(length(selected2[[i]]) != 0){
	 #    		selected2 <- selected2[[i]]
	 #    		break
	 #    	}
	 #    	if(i==1){
	 #    		print(selected2)
	 #    		stop("All selected sets empty")
	 #    	}
	 #    }
	    if(length(selected1) == 0 | length(selected2) == 0){
	    	stop("all(selected1 == 0) | all(selected2 == 0)")
	    }
	    if(length(selected1) == 1 | length(selected2) == 1){
	    	warning("length(selected1) == 1 | length(selected2) == 1)")
	    }
	    # In the (2i - 1)th column of selected_features, make the jth entry
	    # equal 1 if that feature was selected in selected1, and leave it
	    # as zero otherwise. Likewise with the (2i)th column and selected2.
	    selected_features[selected1, 2*i-1] <- 1
	    selected_features[selected2, 2*i] <- 1
	}
	for(i in 1:(2*B)){
		if(all(selected_features[, i] == 0)){
			print(i)
			print(selected_features[, i])
			stop("all(selected_features[, i] == 0)")
		}
	}
	return(selected_features)
}

##### Implementation of complementary pair sampling using apply function to
# consider:

# # Per Samworth & Shah, take a sample (without replacement)
# 		# of size floor(n/2)
# 		samp_size <- floor(n/2)
# 		# p by 2B matrix; in each of 2B columns, the jth row will equal 1
# 		# if that feature was selected, and zero otherwise
# 		selected_features <- matrix(numeric(p*2*B), ncol=(2*B))
# 		# In order to take advantage of apply, create a (floor(n/2)*p +
# 		# floor(n/2) by 2*B matrix. The first floor(n/2)*p elements of each 
# 		# column will be the (unrolled) design matrix, sampled randomly. (Each
# 		# pair of columns will be complementary pairs; there will be B
# 		# complementary pairs yielding 2*B columns). The remaining floor(n/2)
# 		# elements of each columns will be the response y variables. Then this
# 		# matrix will go through an apply function; the relevant function
# 		# will reshape the matrix and fit a lasso to each sample.

# 		apply_matrix <- matrix(numeric(samp_size*(p+1)*2*B), ncol=2*B)
# 		# Take samplings of random orderings
# 		orderings <- matrix(numeric(n*B), ncol=B)
# 		for(i in 1:B){
# 			orderings[, i] <- sample(n)
# 		}
# 		# Form 2*B complementary pair samples of size floor(n/2)
# 		orderings_2B <- matrix(numeric(samp_size*2*B), ncol=2*B)
# 		for(i in 1:B){
# 			orderings_2B[, 2*i - 1] <- orderings[, i][1:samp_size]
# 			orderings_2B[, 2*i] <- orderings[, i][(samp_size + 1):(2*samp_size)]
# 		}
# 		# Fill apply_matrix as described
# 		for(i in 1:(2*B)){
# 			apply_matrix[1:(samp_size*p), i] <- as.numeric(model$
# 				x[orderings_2B[, i], ])
# 			apply_matrix[(samp_size*p+1):nrow(apply_matrix), i] <-
# 				draw[orderings_2B[, i]]
# 		}
# 		# Apply function to apply_matrix and store results in selected_features
# 		results <- apply(apply_matrix, 2, lassoFS, n=samp_size, p=p, q=q)
# 		# results is the transpose of selected_features as I described it.
# 		selected_features <- t(results)


complementaryPairSamplingCoefficients <- function(X.dat, y, B, q){
	# Output: returns a p x 2B matrix of feature coefficients. Each column
	# contains a set of feature coefficients; row i of column j contains the
	# coefficient from the lasso fit on resample j. (If it equals 0, the feature
	# was not selected on resample j.)
	#
	# Inputs:
	# 
	# X.dat: design matrix
	# y: responses
	# B: number of resamples
	# q: Number of features to select on each lasso iteration (glmnet chooses a 
	# suitable lambda such that (approximately) q features are selected)

	n <- nrow(X.dat)
	p <- ncol(X.dat)
	# Per Samworth & Shah, take a sample (without replacement)
	# of size floor(n/2)
	samp_size <- floor(n/2)
	# p by 2*B matrix; in each of B columns, the jth row will equal the 
	# coefficient of the jth feature in the ith reptition.
	feature_coefficients <- matrix(numeric(p*2*B), ncol=(2*B))
	# Complete B repetitions of the following:
	for(i in 1:B){
		# Per Samworth & Shah, take complementary pairwise samples
		order <- sample(n)
		sample1 <- order[1:samp_size]
		sample2 <- order[(samp_size+1):(samp_size+samp_size)]
		# Fit the lasso using glmnet, choosing range of lambdas such that a
		# maximum of q features are selected.
		# use suppressWarnings because otherwise it gives this useless
		# error about q working properly
		fit1 <- suppressWarnings(glmnet(x=X.dat[sample1, ],
			y=y[sample1], family="gaussian", pmax=q))
		fit2 <- suppressWarnings(glmnet(x=X.dat[sample2, ],
			y=y[sample2], family="gaussian", pmax=q))
		# Make the (2*i-1)th column of feature_coefficients equal the
	    # coefficients from lasso fit 1 (choosing the coefficients when q
	    # features are selected), and the (2*i)th column has coefficients
	    # from lasso fit 2.
	    feature_coefficients[, (2*i-1)] <- fit1$beta[, ncol(fit1$beta)]
	    feature_coefficients[, (2*i)] <- fit2$beta[, ncol(fit2$beta)]
	}
	return(feature_coefficients)
}

##### Implementation of complementary pair sampling (returning coefficients)
# using apply function to consider:

# # Per Samworth & Shah, take a sample (without replacement)
# 		# of size floor(n/2)
# 		samp_size <- floor(n/2)
# 		# p by 2*B matrix; in each of B columns, the jth row will equal the 
# 		# coefficient of the jth feature in the ith reptition.
# 		feature_coefficients <- matrix(numeric(p*2*B), ncol=(2*B))
# 		# In order to take advantage of apply, create a (floor(n/2)*p +
# 		# floor(n/2) by 2*B matrix. The first floor(n/2)*p elements of each 
# 		# column will be the (unrolled) design matrix, sampled randomly. (Each
# 		# pair of columns will be complementary pairs; there will be B
# 		# complementary pairs yielding 2*B columns). The remaining floor(n/2)
# 		# elements of each columns will be the response y variables. Then this
# 		# matrix will go through an apply function; the relevant function
# 		# will reshape the matrix and fit a lasso to each sample.

# 		apply_matrix <- matrix(numeric(samp_size*(p+1)*2*B), ncol=2*B)
# 		# Take samplings of random orderings
# 		orderings <- matrix(numeric(n*B), ncol=B)
# 		for(i in 1:B){
# 			orderings[, i] <- sample(n)
# 		}
# 		# Form 2*B complementary pair samples of size floor(n/2)
# 		orderings_2B <- matrix(numeric(samp_size*2*B), ncol=2*B)
# 		for(i in 1:B){
# 			orderings_2B[, 2*i - 1] <- orderings[, i][1:samp_size]
# 			orderings_2B[, 2*i] <- orderings[, i][(samp_size + 1):(2*samp_size)]
# 		}
# 		# Fill apply_matrix as described
# 		for(i in 1:(2*B)){
# 			apply_matrix[1:(samp_size*p), i] <- as.numeric(model$
# 				x[orderings_2B[, i], ])
# 			apply_matrix[(samp_size*p+1):nrow(apply_matrix), i] <-
# 				draw[orderings_2B[, i]]
# 		}
# 		# Apply function to apply_matrix and store results in selected_features
# 		results <- apply(apply_matrix, 2, lassoFSVector, n=samp_size, p=p, q=q)
# 		# results is the transpose of selected_features as I described it.
# 		feature_coefficients <- t(results)

bootstrapSamplingIndicators <- function(X.dat, y, B, q){
	# Outputs: 
	# selected_features: a p x B matrix of selected features. Each column contains
	# a set of selected features, with a 1 in row i of column j if feature i was
	# selected in iteration j, and a 0 otherwise.
	#
	# Inputs:
	# X.dat: design matrix
	# y: responses
	# B: number of resamples
	# q: Number of features to select on each lasso iteration (glmnet chooses a 
	# suitable lambda such that (approximately) q features are selected)

	n <- nrow(X.dat)
	p <- ncol(X.dat)
	# Take a sample (with replacement) of size n
	samp_size <- n
	# p by B matrix; in each of B columns, the jth row will equal 1
	# if that feature was selected, and zero otherwise
	selected_features <- matrix(numeric(p*B), ncol=B)
	# Complete B repetitions of the following:
	for(i in 1:B){
		# Take samples with replacement of size n
		order <- sample(samp_size, size=samp_size, replace=TRUE)
		# Fit the lasso using glmnet, choosing range of lambdas such that a
		# maximum of q features are selected.
		# use suppressWarnings because otherwise it gives this useless
		# error about q working properly
		fit <- suppressWarnings(glmnet(x=X.dat[order, ],
			y=y[order], family="gaussian", pmax=q))
		# fit <- glmnet(x=model$x[order, ],
		# 	y=draw[order], family="gaussian", pmax=q)
		## which coefficients are non-zero?
	    selected <- predict(fit, type="nonzero")
	    # Choose the features chosen at the particular lambda where q
	    # features were selected (i.e. the smallest or last lambda)
	    sels <- selected[[length(selected)]]
	    # In the ith column of selected_features, make the jth entry
	    # equal 1 if that feature was selected in selected, and leave it
	    # as zero otherwise.
	    selected_features[sels, i] <- 1
	}
	return(selected_features)
}

bootstrapSamplingIndicatorsQMax <- function(X.dat, y, B, nlambda){
	# Output: 
	# n_feature_sets: a B-vector that contains the number of selected
	# feature sets for each lasso fit. This ith entry is the number of selected
	# feature sets from the ith lasso fit.
	# selected_features: a p x sum(n_feature_sets) matrix of selected features.
	# Each column contains a set of selected features, with a 1 in row i of
	# column j if feature i was selected in iteration j, and a 0 otherwise.
	# selected_features_lambdas: a sum(n_feature_sets)-vector that
	# keeps track of the lambda corresponding to each selected feature set in
	# selected_features. The ith entry is the lambda corresponding to the
	# selected feature set in the ith column of selected_features.
	#
	# Inputs:
	# X.dat: design matrix
	# y: responses
	# B: number of resamples
	# nlambda: Number of lambdas to use in full sequence
	n <- nrow(X.dat)
	p <- ncol(X.dat)
	# Take a bootstrap sample (with replacement) of size n
	samp_size <- n
	# n_feature_sets is a B-vector that will keep track of the number of selected
	# feature sets for each lasso fit. The ith entry is the number of selected
	# feature sets from the ith lasso fit.
	n_feature_sets <- numeric(B)
	# selected_features will be a p by sum(n_feature_sets) matrix; in each of
	# (sum_i nlambda_i) columns, the jth row will equal 1 if that feature 
	# was selected, and zero otherwise. Get it started with one empty column
	# of the correct size (to be removed at end)
	selected_features <- matrix(0, p, 1)
	# selected_features_lambdas will be a sum(n_feature_sets)-vector that will
	# keep track of the lambda corresponding to each selected feature set in
	# selected_features. The ith entry is the lambda corresponding to the
	# selected feature set in the ith column of selected_features.
	selected_features_lambdas <- numeric()

	# Prepare a samp_size by B matrix of samples with replacement of size n
	# (each of the B columns of the matrix is one sample).
	order <- matrix(sample(n, samp_size*B, replace=TRUE), samp_size, B)
	# Determine lambda sequence to use
	lambda_seq <- sort(glmnet(x=X.dat[order[, 1], ], y=y[order[, 1]],
			family="gaussian", nlambda=nlambda)$lambda, decreasing=TRUE)
	# Complete B repetitions of the following:
	for(i in 1:B){
		# Fit the lasso using glmnet
		fit <- glmnet(x=X.dat[order[, i], ], y=y[order[, i]],
			family="gaussian", lambda=lambda_seq)
		# which coefficients are non-zero?
	    selected <- (as.matrix(fit$beta) != 0)*1
	    # print(selected)
	    # Keep track of the number of selected feature sets for this lasso fit
	    n_feature_sets[i] <- ncol(selected)
	    # print("Number of feature sets selected on this iteration:")
	    # print(n_feature_sets[i])
	    # Combine these feature selections with the previously selected
	    # features
	    selected_features <- cbind(selected_features, selected)
	    # if(length(fit$lambda)==ncol(fit$beta)){
	    	selected_features_lambdas <- c(selected_features_lambdas,
	    		fit$lambda)
    	# }
	    # } 
	    # else{
	    # 	print("ERROR: length(fit$lambda)!=ncol(fit$beta)")
	    # }
	    
	}
	# Remove first empty column of selected features
	selected_features <- selected_features[, 2:ncol(selected_features)]
	# if(!is.matrix(selected_features) || !is.numeric(selected_features)){
	#     print("Error2: selected_features is not a matrix.")
	# }
	if(length(selected_features_lambdas) != ncol(selected_features)){
		print("ERROR: length(selected_features_lambdas) != ncol(selected_features)")
	}
	# print("length(selected_features_lambdas):")
	# print(length(selected_features_lambdas))
	# print("ncol(selected_features):")
	# print(ncol(selected_features))
	return(list(selected_features, n_feature_sets, selected_features_lambdas))
}


bootstrapSamplingCoefficients <- function(X.dat, y, B, q){
	# Output: returns a p x B matrix of feature coefficients. Each column
	# contains a set of feature coefficients; row i of column j contains the
	# coefficient from the lasso fit on resample j. (If it equals 0, the feature
	# was not selected on resample j.)
	#
	# Inputs:
	# X.dat: design matrix
	# y: responses
	# B: number of resamples
	# q: Number of features to select on each lasso iteration (glmnet chooses a 
	# suitable lambda such that (approximately) q features are selected)

	n <- nrow(X.dat)
	p <- ncol(X.dat)
	# Take a sample (with replacement) of size n
	samp_size <- n
	# p by B matrix; in each of B columns, the jth column will contain the 
	# cofficients from the jth lasso fit
	feature_coefficients <- matrix(numeric(p*B), ncol=B)
	# Complete B repetitions of the following:
	for(i in 1:B){
		# Take samples with replacement of size n
		order <- sample(samp_size, size=samp_size, replace=TRUE)
		# Fit the lasso using glmnet, choosing range of lambdas such that a
		# maximum of q features are selected.
		# use suppressWarnings because otherwise it gives this useless
		# error about q working properly
		fit <- suppressWarnings(glmnet(x=X.dat[order, ],
			y=y[order], family="gaussian", pmax=q))
		# Make the ith column of feature_coefficients equal the
	    # coefficients from lasso fit 1 (choosing the coefficients when q
	    # features are selected), and the (2*i)th column has coefficients
	    # from lasso fit 2.
	    feature_coefficients[, i] <- fit$beta[, ncol(fit$beta)]
	}
return(feature_coefficients)
}

####################


# qs used for this iteration

qs <- c(1, 2, 3, 4, 5)

qs_0.51 <- c(1)

### Lasso stability selection--Shah & Samworth technique

lassoSS <- new_method("lassoSS", "Lasso SS (S&S)",
	method = function(model, draw, cutoff, pfer) {
	fit <- stabsel(x=model$x, y=draw, fitfun=glmnet.lasso,
		cutoff=cutoff, q=q, sampling_type="SS")
	return(list(selected=fit$selected))
	},
	# Cutoff probability (proportion of time a feature must be selected by base
	# feature selection method to be included by stability selection) and Per 
	# family error rate (expected number of false selections) for stability
	# selection
	settings = list(cutoff = 0.7, q=10)
)

# make_lasso_ss <- function(q) {
#   	new_method(name = sprintf("lasso_SS%s", q),
#         label = sprintf("SSSS (q=%s)", q),
#         method = function(model, draw, cutoff, q) {
# 		fit <- stabsel(x=model$x, y=draw, fitfun=glmnet.lasso,
# 			cutoff=cutoff, q=q, sampling_type="SS")
# 		return(list(selected=fit$selected))
# 		},
# 		# Cutoff probability (proportion of time a feature must be selected by base
# 		# feature selection method to be included by stability selection) and Per 
# 		# family error rate (expected number of false selections) for stability
# 		# selection
# 		settings = list(cutoff = 0.75, q = q)
# 	)
# }

# # list_of_lasso_sss <- sapply(c(1, 10, 50, 100, 200, 300, 400, 500),
# # 	make_lasso_ss)


make_lasso_ss_columns <- function(mat) {
	# mat is an n x 2 matrix. In each of the n rows, the first term is the
	# the desired cutoff and the second is the desired number of features
	# selected in each lasso fit (q) 
	ret <- list()
	for(i in 1:nrow(mat)){
		meth <- new_method(name = sprintf("lasso_SS%s_%s", mat[i, 1],
			mat[i, 2]),
        label = sprintf("SSSS (tau=%s q=%s)", mat[i, 1], mat[i, 2]),
        method = function(model, draw, cutoff, q) {
			fit <- stabsel(x=model$x, y=draw, fitfun=glmnet.lasso,
				cutoff=cutoff, q=q, sampling_type="SS")
			return(list(selected=fit$selected, beta=NA))
		}, 
		# Cutoff probability (proportion of time a feature must be selected by base
		# feature selection method to be included by stability selection) and q
		# (expected number of features selected on each lasso fit)
		settings = list(cutoff = mat[i, 1], q = mat[i, 2])
		)
		ret <- c(ret, meth)
	}
  	return(ret)
}

lasso_ss_pars_.9 <- matrix(c(rep(0.9, length(qs)), qs), ncol=2, byrow=F)

list_of_lasso_sss_5_.9 <- make_lasso_ss_columns(lasso_ss_pars_.9)

lasso_ss_pars_.7 <- matrix(c(rep(0.7, length(qs)), qs), ncol=2, byrow=F)

list_of_lasso_sss_5_.7 <- make_lasso_ss_columns(lasso_ss_pars_.7)

lasso_ss_pars_.51 <- matrix(c(rep(0.51, length(qs_0.51)), qs_0.51), ncol=2,
	byrow=F)

list_of_lasso_sss_5_.51 <- make_lasso_ss_columns(lasso_ss_pars_.51)


make_lasso_ss_columns_ps <- function(mat) {
	# Returns a p-vector of the estimated probabilities for each feature.
	# 
	# mat is an n x 2 matrix. In each of the n rows, the first term is the
	# the desired cutoff and the second is the desired number of features
	# selected in each lasso fit (q) 
	ret <- list()
	for(i in 1:nrow(mat)){
		meth <- new_method(name = sprintf("lasso_SS%s_%s", mat[i, 1],
			mat[i, 2]),
        label = sprintf("SSSS (tau=%s q=%s)", mat[i, 1], mat[i, 2]),
        method = function(model, draw, cutoff, q) {
			fit <- stabsel(x=model$x, y=draw, fitfun=glmnet.lasso,
				cutoff=cutoff, q=q, sampling_type="SS")
			return(list(phat=fit$phat[, ncol(fit$phat)], beta=NA))
		}, 
		# Cutoff probability (proportion of time a feature must be selected by base
		# feature selection method to be included by stability selection) and q
		# (expected number of features selected on each lasso fit)
		settings = list(cutoff = mat[i, 1], q = mat[i, 2])
		)
		ret <- c(ret, meth)
	}
  	return(ret)
}

lasso_ss_ps_pars_.7 <- matrix(c(rep(0.7, length(qs)), qs), ncol=2, byrow=F)

list_of_lasso_ps_sss_5_.7 <- make_lasso_ss_columns_ps(lasso_ss_ps_pars_.7)

lassoSS_phat <- new_method("lassoSS_phat", "Lasso SS (S&S) p_hats",
	method = function(model, draw, cutoff, q, pfer, B) {
		R <- diag(ncol(model$x))
	fit <- stabselGreg(x=model$x, y=draw, fitfun=glmnet.lasso,
		cutoff=cutoff, q=q, B=B, sampling_type="SS"
		, R=R
		)
	# return(list(phat=fit$phat[, ncol(fit$phat)], beta=NA))
	return(list(phat=fit$feat_sel_props, beta=NA,
		first_lasso_selected=fit$feat_sel_mat[1, ]))
	},
	# Cutoff probability (proportion of time a feature must be selected by base
	# feature selection method to be included by stability selection) and Per 
	# family error rate (expected number of false selections) for stability
	# selection
	settings = list(cutoff = 0.7, q=11, B=100)
)

SS_SS_random <- new_method("SS_SS_random",
	"S&S SS (random X)",
	method = function(model, draw, cutoff, pfer, B) {
	require(glmnet)
	R <- diag(ncol(draw$X))
	n <- nrow(draw$X)
	if(length(draw$y) != n){
		stop("length(draw$y) != n")
	}
	# Determine q: do lasso with cross-validation on full data sample.
	# Cross-validated model size will be the size we use for stability
	# selection.
	inds_size <- sample(1:n, floor(n/2))
	X <- draw$X[inds_size, ]
	y <- draw$y[inds_size]
	colnames(X) <- character()
	rownames(X) <- character()
	names(y) <- character()
	size_results <- cv.glmnet(x=X, y=y, parallel=TRUE, family="gaussian")
	rm(X)
	rm(y)
	q <- nrow(predict(size_results, s = "lambda.min", type="nonzero"))
	fit <- stabselGreg(x=draw$X, y=draw$y, fitfun=glmnet.lasso,
		cutoff=cutoff, q=q, B=B, sampling_type="SS"
		, R=R
		)
	# return(list(phat=fit$phat[, ncol(fit$phat)], beta=NA))
	return(list(phat=fit$feat_sel_props, beta=NA,
		first_lasso_selected=fit$feat_sel_mat[1, ]))
	},
	# Cutoff probability (proportion of time a feature must be selected by base
	# feature selection method to be included by stability selection) and Per 
	# family error rate (expected number of false selections) for stability
	# selection
	settings = list(cutoff = 0.7, B=100)
)

SS_SS_cssr <- new_method("SS_SS_cssr",
	"S&S SS (random X, custom SS function)",
	method = function(model, draw, B, model_size) {
	# Original stability selection procedure
	# Get lambda
	lambda <- cssr::getLassoLambda(X=draw$X, y=draw$y, lambda_choice="min")

	# Don't provide clusters
	res <- cssr::css(X=draw$X, y=draw$y, lambda=lambda,
		num_cores=detectCores() - 1)

	# Confirm no clusters in the results
	stopifnot(ncol(res$clus_sel_mat) == ncol(draw$X))

	selected <- list()
	selected_clusts <- list()
	for(i in 1:model_size){
		res_i <- cssr::getCssSelections(res, min_num_clusts=i,
			max_num_clusts=i)

		set_i <- res_i$selected_feats
		clusts_i <- res_i$selected_clusts

		# Number of clusters should equal number of features for regular
		# stability selection (each "cluster" should contain a single feature)
		stopifnot(length(set_i) == length(clusts_i))

		if(length(set_i) == i){
			selected[[i]] <- set_i
			selected_clusts[[i]] <- clusts_i
		}
	}
	# weighting doesn't matter since no clusters provided--using any weighting
	# scheme will yield the same results in getCssPreds or getCssSelections
	return(list(css_res=res, selected=selected, selected_clusts=selected_clusts,
		method="sparse", testX=draw$testX, testY=draw$testY,
		testMu=draw$testMu))
	},
	settings = list(B=100, model_size=11)
)

SS_SS_random_custom <- new_method("SS_SS_random_custom",
	"S&S SS (random X, custom SS function)",
	method = function(model, draw, B) {
	# Original stability selection procedure
	require(glmnet)
	R <- diag(ncol(draw$X))
	n <- nrow(draw$X)
	if(length(draw$y) != n){
		stop("length(draw$y) != n")
	}
	# Determine q: do lasso with cross-validation on full data sample.
	# Cross-validated model size will be the size we use for stability
	# selection.
	inds_size <- sample(1:n, floor(n/2))
	X <- draw$X[inds_size, ]
	y <- draw$y[inds_size]
	colnames(X) <- character()
	rownames(X) <- character()
	names(y) <- character()
	size_results <- cv.glmnet(x=X, y=y, parallel=TRUE, family="gaussian")
	rm(X)
	rm(y)

	# Get selection proportions
	css_results <- css(X=draw$X, y=draw$y, lambda=size_results$lambda.min,
		B=B, sampling_type="SS")

	phat <- colMeans(css_results$feat_sel_mat)

	# Get clusters and sort in decreasing order of selection proportion
	clus_sel_props <- colMeans(css_results$clus_sel_mat)
	sel_clusts <- css_results$clusters[names(sort(clus_sel_props, decreasing=TRUE))]

	selected <- getSelectionPrototypes(css_results, sel_clusts)

	# fit <- cssWrapper(X=draw$X, y=draw$y, lambda=size_results$lambda.min,
	# 	B=B, sampling_type="SS")


	# return(list(phat=fit$phat[, ncol(fit$phat)], beta=NA))
	return(list(phat=phat, beta=NA,
		first_lasso_selected=css_results$feat_sel_mat[1, ], selected=selected,
		tied_with_next=getTiedWithNext(names(sel_clusts), clus_sel_props)
		))
	# return(list(phat=fit$feat_sel_props, beta=NA,
	# 	first_lasso_selected=fit$feat_sel_mat[1, ], selected=fit$selected,
	# 	tied_with_next=getTiedWithNext(names(fit$selected_clusts),
	# 		fit$clus_sel_props)
	# 	))
	},
	# Cutoff probability (proportion of time a feature must be selected by base
	# feature selection method to be included by stability selection) and Per 
	# family error rate (expected number of false selections) for stability
	# selection
	settings = list(B=100)
)

lassoMB_phat <- new_method("lassoMB_phat", "Lasso SS (M&B) p_hats",
	method = function(model, draw, cutoff, q, pfer, B) {
		R <- diag(ncol(model$x))
	fit <- stabselGreg(x=model$x, y=draw, fitfun=glmnet.lasso,
		cutoff=cutoff, q=q, B=B, sampling_type="MB"
		, R=R
		)
	# return(list(phat=fit$phat[, ncol(fit$phat)], beta=NA))
	return(list(phat=fit$feat_sel_props, beta=NA,
		first_lasso_selected=fit$feat_sel_mat[1, ]))
	},
	# Cutoff probability (proportion of time a feature must be selected by base
	# feature selection method to be included by stability selection) and Per 
	# family error rate (expected number of false selections) for stability
	# selection
	settings = list(cutoff = 0.7, q=11, B=100)
)

lassoSS_phat_ideal <- new_method("lassoSS_phat_ideal",
	"Ideal Lasso SS (S&S) p_hats",
	method = function(model, draw, cutoff, q, pfer, B) {
		# Set R to exactly reflect actual (known) clusters
		R <- (model$Sigma !=0)*1
		# Check if used a model where signal cluster variable was removed; if
		# so, remove first row and column of Sigma (and therefore R)
		# print(str(R))
		# print(str(model$x))
		# if(ncol(R) > ncol(model$x)){
		# 	R <- R[2:nrow(R), 2:ncol(R)]
		# }
	fit <- stabselGreg(x=model$x, y=draw, fitfun=glmnet.lasso,
		cutoff=cutoff, q=q, B=B, sampling_type="SS"
		, R=R
		)
	# return(list(phat=fit$phat[, ncol(fit$phat)], beta=NA))
	return(list(phat=fit$feat_sel_props, beta=NA, R=R, selected=fit$selected))
	},
	# Cutoff probability (proportion of time a feature must be selected by base
	# feature selection method to be included by stability selection) and Per 
	# family error rate (expected number of false selections) for stability
	# selection
	settings = list(cutoff = 0.52, q=11, B=100)
)

SS_GSS_random <- new_method("SS_GSS_random",
	"S&S GSS (random X)",
	method = function(model, draw, cutoff, pfer, B) {
	require(glmnet)
	# Take out first row of Sigma that no longer needs to be there (again, 
    # assuming one block)
    if(model$nblocks==1 & model$sig_blocks ==1){
        Sigma <- model$Sigma[2:nrow(model$Sigma), 2:ncol(model$Sigma)]
    }
	# Set R to exactly reflect actual (known) clusters
	R <- (Sigma !=0)*1

	n <- nrow(draw$X)
	if(length(draw$y) != n){
		stop("length(draw$y) != n")
	}

	
	# Determine q: do lasso with cross-validation on full data sample.
	# Cross-validated model size will be the size we use for stability
	# selection.
	inds_size <- sample(1:n, floor(n/2))
	size_results <- cv.glmnet(x=draw$X[inds_size, ], y=draw$y[inds_size],
		parallel=TRUE, family="gaussian")
	q <- nrow(predict(size_results, s = "lambda.min", type="nonzero"))

	# Check if used a model where signal cluster variable was removed; if
	# so, remove first row and column of Sigma (and therefore R)
	# print(str(R))
	# print(str(model$x))
	# if(ncol(R) > ncol(model$x)){
	# 	R <- R[2:nrow(R), 2:ncol(R)]
	# }
	fit <- stabselGreg(x=draw$X, y=draw$y, fitfun=glmnet.lasso,
		cutoff=cutoff, q=q, B=B, sampling_type="SS"
		, R=R
		)
	# return(list(phat=fit$phat[, ncol(fit$phat)], beta=NA))
	return(list(phat=fit$feat_sel_props, beta=NA, R=R, selected=fit$selected))
	},
	# Cutoff probability (proportion of time a feature must be selected by base
	# feature selection method to be included by stability selection) and Per 
	# family error rate (expected number of false selections) for stability
	# selection
	settings = list(cutoff = 0.52, B=100)
)

SS_CSS_sparse_cssr <- new_method("SS_CSS_sparse_cssr",
	"S&S GSS (random X, custom SS function)",
	method = function(model, draw, B, model_size) {
	# Sparse cluster stability selection
	lambda <- cssr::getLassoLambda(draw$X, draw$y, lambda_choice="min")

	stopifnot(model$nblocks == 1)
	stopifnot(model$sig_blocks == 1)

	res <- cssr::css(draw$X, draw$y, lambda=lambda, clusters=1:model$block_size,
		num_cores=detectCores() - 1)

	# Confirm cluster in the results
	stopifnot(ncol(res$clus_sel_mat) == ncol(draw$X) - model$block_size + 1)

	selected <- list()
	selected_clusts <- list()
	for(i in 1:model_size){
		res_i <- cssr::getCssSelections(res, weighting="sparse",
			min_num_clusts=i, max_num_clusts=i)

		set_i <- res_i$selected_feats
		clusts_i <- res_i$selected_clusts

		stopifnot(length(clusts_i) <= length(set_i))
		
		if(length(clusts_i) == i){
			selected[[i]] <- set_i
			selected_clusts[[i]] <- clusts_i
		}
	}

	return(list(css_res=res, selected=selected, selected_clusts=selected_clusts,
		method="sparse", testX=draw$testX, testY=draw$testY,
		testMu=draw$testMu))

	},
	settings = list(B=100, model_size=11)
)




SS_GSS_random_custom <- new_method("SS_GSS_random_custom",
	"S&S GSS (random X, custom SS function)",
	method = function(model, draw, B) {
	require(glmnet)
	# Sparse cluster stability selection

	# Take out first row of Sigma that no longer needs to be there (again, 
    # assuming one block)
    if(model$nblocks==1 & model$sig_blocks ==1){
        Sigma <- model$Sigma[2:nrow(model$Sigma), 2:ncol(model$Sigma)]
    }
	# Set R to exactly reflect actual (known) clusters
	R <- (Sigma !=0)*1

	n <- nrow(draw$X)
	if(length(draw$y) != n){
		stop("length(draw$y) != n")
	}

	
	# Determine q: do lasso with cross-validation on full data sample.
	# Cross-validated model size will be the size we use for stability
	# selection.
	inds_size <- sample(1:n, floor(n/2))
	size_results <- cv.glmnet(x=draw$X[inds_size, ], y=draw$y[inds_size],
		parallel=TRUE, family="gaussian")

	clusters <- getClustersFromR(R)


	# Get selection proportions
	css_results <- css(X=draw$X, y=draw$y, lambda=size_results$lambda.min,
		clusters=clusters, B=B, sampling_type="SS")

	phat <- colMeans(css_results$feat_sel_mat)

	# Get clusters and sort in decreasing order of selection proportion
	clus_sel_props <- colMeans(css_results$clus_sel_mat)
	sel_clusts <- css_results$clusters[names(sort(clus_sel_props,
		decreasing=TRUE))]

	selected <- getSelectionPrototypes(css_results, sel_clusts)

	# fit <- cssWrapper(X=draw$X, y=draw$y, lambda=size_results$lambda.min,
	# 	B=B, sampling_type="SS")


	# return(list(phat=fit$phat[, ncol(fit$phat)], beta=NA))
	return(list(phat=phat, beta=NA,
		first_lasso_selected=css_results$feat_sel_mat[1, ], selected=selected,
		tied_with_next=getTiedWithNext(names(sel_clusts), clus_sel_props)
		))





	# fit <- cssWrapper(X=draw$X, y=draw$y, lambda=size_results$lambda.min, B=B,
	# 	clusters=clusters, sampling_type="SS", weighting="sparse")
	# # return(list(phat=fit$phat[, ncol(fit$phat)], beta=NA))
	# return(list(phat=fit$feat_sel_props, beta=NA, R=R, selected=fit$selected,
	# 	tied_with_next=getTiedWithNext(names(fit$selected_clusts),
	# 		fit$clus_sel_props)
	# 	))
	},
	# Cutoff probability (proportion of time a feature must be selected by base
	# feature selection method to be included by stability selection) and Per 
	# family error rate (expected number of false selections) for stability
	# selection
	settings = list(B=100)
)

SS_GSS_random_avg <- new_method("SS_GSS_random_avg",
	"S&S GSS (random X, averaging)",
	method = function(model, draw, cutoff, pfer, B) {
	require(glmnet)
	# Outputs:
    #
    # to_avg: A logical vector to_avg of the same length as 
    # selected. If feature j is in a cluster, the jth entry of to_avg will be
    # TRUE. 
    #
    # selected_clusts: A list of the same length as selected. If feature j is in a
    # cluster, the jth entry of selected_clusts will contain the features from the
    # cluster of which feature j is a member.
    #
    # weights: A list of the same length as selected. If feature j is in a
    # cluster, the jth entry of weights will
    # be the weights to use (in the same order as the jth entry of
    # selected_clusts).

	# Take out first row of Sigma that no longer needs to be there (again, 
    # assuming one block)
    if(model$nblocks==1 & model$sig_blocks ==1){
        Sigma <- model$Sigma[2:nrow(model$Sigma), 2:ncol(model$Sigma)]
    }
	# Set R to exactly reflect actual (known) clusters
	R <- (Sigma !=0)*1

	n <- nrow(draw$X)
	if(length(draw$y) != n){
		stop("length(draw$y) != n")
	}
	# Determine q: do lasso with cross-validation on full data sample.
	# Cross-validated model size will be the size we use for stability
	# selection.
	inds_size <- sample(1:n, floor(n/2))
	size_results <- cv.glmnet(x=draw$X[inds_size, ], y=draw$y[inds_size],
		parallel=TRUE, family="gaussian")
	q <- nrow(predict(size_results, s = "lambda.min", type="nonzero"))
	# Check if used a model where signal cluster variable was removed; if
	# so, remove first row and column of Sigma (and therefore R)
	# print(str(R))
	# print(str(model$x))
	# if(ncol(R) > ncol(model$x)){
	# 	R <- R[2:nrow(R), 2:ncol(R)]
	# }
	fit <- stabselGreg(x=draw$X, y=draw$y, fitfun=glmnet.lasso,
		cutoff=cutoff, q=q, B=B, sampling_type="SS"
		, R=R, average=TRUE
		)

	# print("to_avg:")
	# print(fit$to_avg)
	# print("selected_clusts:")
	# print(fit$selected_clusts)
	# print("weights:")
	# print(fit$weights)
	# return(list(phat=fit$phat[, ncol(fit$phat)], beta=NA))
	return(list(phat=fit$feat_sel_props, beta=NA, R=R, selected=fit$selected,
		to_avg = fit$to_avg, selected_clusts = fit$selected_clusts, weights = fit$weights))
	},
	# Cutoff probability (proportion of time a feature must be selected by base
	# feature selection method to be included by stability selection) and Per 
	# family error rate (expected number of false selections) for stability
	# selection
	settings = list(cutoff = 0.52, B=100)
)


SS_CSS_weighted_cssr <- new_method("SS_CSS_weighted_cssr",
	"S&S GSS (random X, averaging, custom SS function)",
	method = function(model, draw, B, model_size) {
	# Sparse cluster stability selection
	lambda <- cssr::getLassoLambda(draw$X, draw$y, lambda_choice="min")

	stopifnot(model$nblocks == 1)
	stopifnot(model$sig_blocks == 1)

	res <- cssr::css(draw$X, draw$y, lambda=lambda, clusters=1:model$block_size,
		num_cores=detectCores() - 1)

	# Confirm cluster in the results
	stopifnot(ncol(res$clus_sel_mat) == ncol(draw$X) - model$block_size + 1)

	selected <- list()
	selected_clusts <- list()
	for(i in 1:model_size){
		res_i <- cssr::getCssSelections(res, weighting="weighted_avg",
			min_num_clusts=i, max_num_clusts=i)

		set_i <- res_i$selected_feats
		clusts_i <- res_i$selected_clusts

		stopifnot(length(clusts_i) <= length(set_i))
		
		if(length(clusts_i) == i){
			selected[[i]] <- set_i
			selected_clusts[[i]] <- clusts_i
		}
	}

	return(list(css_res=res, selected=selected, selected_clusts=selected_clusts,
		method="weighted_avg", testX=draw$testX, testY=draw$testY,
		testMu=draw$testMu))

	},
	settings = list(B=100, model_size=11)
)


SS_GSS_random_avg_custom <- new_method("SS_GSS_random_avg_custom",
	"S&S GSS (random X, averaging, custom SS function)",
	method = function(model, draw, B) {
	# Weighted averaged cluster stability selection
	require(glmnet)
	# Outputs:
    #
    # to_avg: A logical vector to_avg of the same length as 
    # selected. If feature j is in a cluster, the jth entry of to_avg will be
    # TRUE. 
    #
    # selected_clusts: A list of the same length as selected. If feature j is in a
    # cluster, the jth entry of selected_clusts will contain the features from the
    # cluster of which feature j is a member.
    #
    # weights: A list of the same length as selected. If feature j is in a
    # cluster, the jth entry of weights will
    # be the weights to use (in the same order as the jth entry of
    # selected_clusts).

	# Take out first row of Sigma that no longer needs to be there (again, 
    # assuming one block)
    if(model$nblocks==1 & model$sig_blocks ==1){
        Sigma <- model$Sigma[2:nrow(model$Sigma), 2:ncol(model$Sigma)]
    }
	# Set R to exactly reflect actual (known) clusters
	R <- (Sigma !=0)*1

	n <- nrow(draw$X)
	if(length(draw$y) != n){
		stop("length(draw$y) != n")
	}
	# Determine q: do lasso with cross-validation on full data sample.
	# Cross-validated model size will be the size we use for stability
	# selection.
	inds_size <- sample(1:n, floor(n/2))
	size_results <- cv.glmnet(x=draw$X[inds_size, ], y=draw$y[inds_size],
		parallel=TRUE, family="gaussian")

	clusters <- getClustersFromR(R)

	# Get selection proportions
	css_results <- css(X=draw$X, y=draw$y, lambda=size_results$lambda.min,
		clusters=clusters, B=B, sampling_type="SS")

	phat <- colMeans(css_results$feat_sel_mat)

	# Get clusters and sort in decreasing order of selection proportion
	clus_sel_props <- colMeans(css_results$clus_sel_mat)
	sorted_clus_sel_props <- sort(clus_sel_props, decreasing=TRUE)
	sel_clusts <- css_results$clusters[names(sorted_clus_sel_props)]

	selected <- getSelectionPrototypes(css_results, sel_clusts)

	# fit <- cssWrapper(X=draw$X, y=draw$y, lambda=size_results$lambda.min,
	# 	B=B, sampling_type="SS")


	# return(list(phat=fit$phat[, ncol(fit$phat)], beta=NA))
	return(list(phat=phat, beta=NA,
		first_lasso_selected=css_results$feat_sel_mat[1, ], selected=selected,
		tied_with_next=getTiedWithNext(names(sel_clusts), clus_sel_props),
		selected_clusts = sel_clusts,
		weights = clustWeights(css_results, sorted_clus_sel_props,
			weighting="weighted_avg")
		))







	# fit <- cssWrapper(X=draw$X, y=draw$y, lambda=size_results$lambda.min, B=B,
	# 	clusters=clusters, sampling_type="SS", weighting="weighted_avg")

 #    stopifnot(is.list(fit$selected_clusts))
 #    stopifnot(length(fit$selected_clusts) <= length(fit$selected))

 #    stopifnot(is.list(fit$weights))
 #    stopifnot(length(fit$selected_clusts) == length(fit$weights))

 #    if(length(fit$selected_clusts) > 0){
 #        for(i in 1:length(fit$selected_clusts)){
 #            stopifnot(is.integer(fit$selected_clusts[[i]]))
 #            stopifnot(length(fit$selected_clusts[[i]]) >= 1)
 #            stopifnot(all(!is.na(fit$selected_clusts[[i]])))
 #            stopifnot(length(unique(fit$selected_clusts[[i]])) ==
 #                length(fit$selected_clusts[[i]]))
 #            stopifnot(fit$selected[i] %in% fit$selected_clusts[[i]])

 #            stopifnot(length(fit$selected_clusts[[i]]) == length(fit$weights[[i]]))
 #            stopifnot(is.numeric(fit$weights[[i]]))
 #            stopifnot(all(fit$weights[[i]] >= 0))
 #            stopifnot(all(fit$weights[[i]] <= 1))
 #            stopifnot(abs(sum(fit$weights[[i]]) - 1) < 10^(-6))
 #        }
 #    }
	
	# return(list(phat=fit$feat_sel_props, beta=NA, R=R, selected=fit$selected,
	# 	to_avg = rep(TRUE, length(fit$selected)),
	# 	selected_clusts = fit$selected_clusts,
	# 	weights = fit$weights,
	# 	tied_with_next=getTiedWithNext(names(fit$selected_clusts),
	# 		fit$clus_sel_props)
	# 	))
	},
	# Cutoff probability (proportion of time a feature must be selected by base
	# feature selection method to be included by stability selection) and Per 
	# family error rate (expected number of false selections) for stability
	# selection
	settings = list(B=100)
)



SS_GSS_random_avg_unwt <- new_method("SS_GSS_random_avg_unwt",
	"S&S GSS (random X, unweighted averaging)",
	method = function(model, draw, cutoff, pfer, B) {
	require(glmnet)
	# Outputs:
    #
    # to_avg: A logical vector to_avg of the same length as 
    # selected. If feature j is in a cluster, the jth entry of to_avg will be
    # TRUE. 
    #
    # selected_clusts: A list of the same length as selected. If feature j is in a
    # cluster, the jth entry of selected_clusts will contain the features from the
    # cluster of which feature j is a member.
    #
    # weights: A list of the same length as selected. If feature j is in a
    # cluster, the jth entry of weights will
    # be the weights to use (in the same order as the jth entry of
    # selected_clusts).

	# Take out first row of Sigma that no longer needs to be there (again, 
    # assuming one block)
    if(model$nblocks==1 & model$sig_blocks ==1){
        Sigma <- model$Sigma[2:nrow(model$Sigma), 2:ncol(model$Sigma)]
    }
	# Set R to exactly reflect actual (known) clusters
	R <- (Sigma !=0)*1

	n <- nrow(draw$X)
	if(length(draw$y) != n){
		stop("length(draw$y) != n")
	}
	# Determine q: do lasso with cross-validation on full data sample.
	# Cross-validated model size will be the size we use for stability
	# selection.
	inds_size <- sample(1:n, floor(n/2))
	size_results <- cv.glmnet(x=draw$X[inds_size, ], y=draw$y[inds_size],
		parallel=TRUE, family="gaussian")
	q <- nrow(predict(size_results, s = "lambda.min", type="nonzero"))
	# Check if used a model where signal cluster variable was removed; if
	# so, remove first row and column of Sigma (and therefore R)
	# print(str(R))
	# print(str(model$x))
	# if(ncol(R) > ncol(model$x)){
	# 	R <- R[2:nrow(R), 2:ncol(R)]
	# }
	fit <- stabselGreg(x=draw$X, y=draw$y, fitfun=glmnet.lasso,
		cutoff=cutoff, q=q, B=B, sampling_type="SS"
		, R=R, average=TRUE, weighted=FALSE
		)

	# print("to_avg:")
	# print(fit$to_avg)
	# print("selected_clusts:")
	# print(fit$selected_clusts)
	# print("weights:")
	# print(fit$weights)
	# return(list(phat=fit$phat[, ncol(fit$phat)], beta=NA))
	return(list(phat=fit$feat_sel_props, beta=NA, R=R, selected=fit$selected,
		to_avg = fit$to_avg, selected_clusts = fit$selected_clusts, weights = fit$weights))
	},
	# Cutoff probability (proportion of time a feature must be selected by base
	# feature selection method to be included by stability selection) and Per 
	# family error rate (expected number of false selections) for stability
	# selection
	settings = list(cutoff = 0.52, B=100)
)

SS_CSS_avg_cssr <- new_method("SS_CSS_avg_cssr",
	"S&S GSS (random X, unweighted averaging, custom SS function)",
	method = function(model, draw, B, model_size) {
	# Simple averaged cluster stability selection
	lambda <- cssr::getLassoLambda(draw$X, draw$y, lambda_choice="min")

	stopifnot(model$nblocks == 1)
	stopifnot(model$sig_blocks == 1)

	res <- cssr::css(draw$X, draw$y, lambda=lambda, clusters=1:model$block_size,
		num_cores=detectCores() - 1)

	# Confirm cluster in the results
	stopifnot(ncol(res$clus_sel_mat) == ncol(draw$X) - model$block_size + 1)

	selected <- list()
	selected_clusts <- list()
	for(i in 1:model_size){
		res_i <- cssr::getCssSelections(res, weighting="simple_avg",
			min_num_clusts=i, max_num_clusts=i)

		set_i <- res_i$selected_feats
		clusts_i <- res_i$selected_clusts

		stopifnot(length(clusts_i) <= length(set_i))
		
		if(length(clusts_i) == i){
			selected[[i]] <- set_i
			selected_clusts[[i]] <- clusts_i
		}
	}

	return(list(css_res=res, selected=selected, selected_clusts=selected_clusts,
		method="simple_avg", testX=draw$testX, testY=draw$testY,
		testMu=draw$testMu))

	},
	settings = list(B=100, model_size=11)
)




SS_GSS_random_avg_unwt_custom <- new_method("SS_GSS_random_avg_unwt_custom",
	"S&S GSS (random X, unweighted averaging, custom SS function)",
	method = function(model, draw, B) {
	# Simple averaged cluster stability selection
	require(glmnet)
	# Outputs:
    #
    # to_avg: A logical vector to_avg of the same length as 
    # selected. If feature j is in a cluster, the jth entry of to_avg will be
    # TRUE. 
    #
    # selected_clusts: A list of the same length as selected. If feature j is in a
    # cluster, the jth entry of selected_clusts will contain the features from the
    # cluster of which feature j is a member.
    #
    # weights: A list of the same length as selected. If feature j is in a
    # cluster, the jth entry of weights will
    # be the weights to use (in the same order as the jth entry of
    # selected_clusts).

	# Take out first row of Sigma that no longer needs to be there (again, 
    # assuming one block)
    if(model$nblocks==1 & model$sig_blocks ==1){
        Sigma <- model$Sigma[2:nrow(model$Sigma), 2:ncol(model$Sigma)]
    }
	# Set R to exactly reflect actual (known) clusters
	R <- (Sigma !=0)*1

	n <- nrow(draw$X)
	if(length(draw$y) != n){
		stop("length(draw$y) != n")
	}
	# Determine q: do lasso with cross-validation on full data sample.
	# Cross-validated model size will be the size we use for stability
	# selection.
	inds_size <- sample(1:n, floor(n/2))
	size_results <- cv.glmnet(x=draw$X[inds_size, ], y=draw$y[inds_size],
		parallel=TRUE, family="gaussian")

	clusters <- getClustersFromR(R)

	clusters <- getClustersFromR(R)

	# Get selection proportions
	css_results <- css(X=draw$X, y=draw$y, lambda=size_results$lambda.min,
		clusters=clusters, B=B, sampling_type="SS")

	phat <- colMeans(css_results$feat_sel_mat)

	# Get clusters and sort in decreasing order of selection proportion
	clus_sel_props <- colMeans(css_results$clus_sel_mat)
	sorted_clus_sel_props <- sort(clus_sel_props, decreasing=TRUE)
	sel_clusts <- css_results$clusters[names(sorted_clus_sel_props)]

	selected <- getSelectionPrototypes(css_results, sel_clusts)

	# return(list(phat=fit$phat[, ncol(fit$phat)], beta=NA))
	return(list(phat=phat, beta=NA,
		first_lasso_selected=css_results$feat_sel_mat[1, ], selected=selected,
		tied_with_next=getTiedWithNext(names(sel_clusts), clus_sel_props),
		selected_clusts = sel_clusts,
		weights = clustWeights(css_results, sorted_clus_sel_props,
			weighting="simple_avg")
		))






	# fit <- cssWrapper(X=draw$X, y=draw$y, lambda=size_results$lambda.min, B=B,
	# 	clusters=clusters, sampling_type="SS", weighting="simple_avg")

 #    stopifnot(is.list(fit$selected_clusts))
 #    stopifnot(length(fit$selected_clusts) <= length(fit$selected))

 #    stopifnot(is.list(fit$weights))
 #    stopifnot(length(fit$selected_clusts) == length(fit$weights))

 #    if(length(fit$selected_clusts) > 0){
 #        for(i in 1:length(fit$selected_clusts)){
 #            stopifnot(is.integer(fit$selected_clusts[[i]]))
 #            stopifnot(length(fit$selected_clusts[[i]]) >= 1)
 #            stopifnot(all(!is.na(fit$selected_clusts[[i]])))
 #            stopifnot(length(unique(fit$selected_clusts[[i]])) ==
 #                length(fit$selected_clusts[[i]]))
 #            stopifnot(fit$selected[i] %in% fit$selected_clusts[[i]])

 #            stopifnot(length(fit$selected_clusts[[i]]) == length(fit$weights[[i]]))
 #            stopifnot(is.numeric(fit$weights[[i]]))
 #            stopifnot(all(fit$weights[[i]] >= 0))
 #            stopifnot(all(fit$weights[[i]] <= 1))
 #            stopifnot(abs(sum(fit$weights[[i]]) - 1) < 10^(-6))

 #        }
 #    }
	
	# return(list(phat=fit$feat_sel_props, beta=NA, R=R, selected=fit$selected,
	# 	to_avg = rep(TRUE, length(fit$selected)),
	# 	selected_clusts = fit$selected_clusts, weights = fit$weights,
	# 	tied_with_next=getTiedWithNext(names(fit$selected_clusts),
	# 		fit$clus_sel_props)
	# 	))
	},
	settings = list(B=100)
)


clusRepLasso_cssr <- new_method("clusRepLasso_cssr",
	"BRVZ method (random X, unweighted averaging)",
	method = function(model, draw, model_size) {
	
	stopifnot(model$nblocks == 1)
	stopifnot(model$sig_blocks == 1)

	res <- cssr::clusterRepLasso(draw$X, draw$y, clusters=1:model$block_size,
		nlambda=2000)

	num_sets <- min(length(res$selected_clusts_list), model_size)

	return(list(selected_clusts_list=res$selected_clusts_list[1:num_sets],
		testX=draw$testX, testY=draw$testY, testMu=draw$testMu))	
	},
	settings = list(model_size=11)
)




BRVZ_avg_unwt <- new_method("BRVZ_avg_unwt",
	"BRVZ method (random X, unweighted averaging)",
	method = function(model, draw) {
	require(glmnet)
	# Outputs:
    #
    # selected_sets: A list of integer vectors. The jth entry is the selected
    # set of size j (counting choosing a cluster average as one feature).
    #
    # to_avg_list: A list of logical vectors. The jth entrry is aa logical vector 
    # of length j. If feature k is in a cluster, the kth entry of the jth vector
    # in to_avg will be TRUE. 
    #
    # selected_clusts_list: A list of lists. The jth entry is a list of length j.
    # If feature k is in a cluster, the kth entry of the jth list in selected_clusts
    # will contain the features from the cluster of which feature k is a member.
    #
    # weights_list: A list of lists. The jth entry is a list of length j.
    # If feature k is in a cluster, the kth entry of the jth list in
    # weights_list will be the weights to use (in the same order as the kth 
    # entry of selected_clusts).

	# Take out first row of Sigma that no longer needs to be there (again, 
    # assuming one block)
    if(model$nblocks==1 & model$sig_blocks ==1){
        Sigma <- model$Sigma[2:nrow(model$Sigma), 2:ncol(model$Sigma)]
    } else{
    	stop("!(model$nblocks==1 & model$sig_blocks ==1)")
    }
	# Set R to exactly reflect actual (known) clusters
	R <- (Sigma !=0)*1

	# Construct design matrix: replace cluster with average of features,
	# leave other features alone

	# # Identify clustered features: row i in R contains all the features
 #    # that are in a cluster with feature i, so come up with a list containing
 #    # all the clusters and then remove the repeats
 #    # print("R:")
 #    # print(R)
 #    clusters <- list()
 #    for(i in 1:nrow(R)){
 #        clusters[[i]] <- which(R[i, ] > 0)
 #    }

 #    clusters <- unique(clusters)
 #    # Only care about clusters with more than one element (only ones that need
 #    # to be treated differently)
 #    # keep track of whether there's more than one cluster or not
 #    multiple <- sum(lengths(clusters) > 1) > 1

 	BRVZ_results <- clusterRepLasso(x=draw$X, y=draw$y, R=R)

	return(list(selected_sets=BRVZ_results$selected_sets,
		to_avg_list=BRVZ_results$to_avg_list,
		selected_clusts_list=BRVZ_results$selected_clusts_list,
		weights_list=BRVZ_results$weights_list))
	}
)

protolasso_cssr <- new_method("protolasso_cssr",
	"Lasso (cluster prototype, Random design)", method=function(model, draw,
		model_size){

	stopifnot(model$nblocks == 1)
	stopifnot(model$sig_blocks == 1)

	res <- cssr::protolasso(draw$X, draw$y, clusters=1:model$block_size,
		nlambda=2000)

	num_sets <- min(length(res$selected_sets), model_size)

	return(list(selected_sets=res$selected_sets[1:num_sets],
		testX=draw$testX, testY=draw$testY, testMu=draw$testMu))	
	},
	settings = list(model_size=11)
)



lasso_proto <- new_method("lasso_proto",
	"Lasso (cluster prototype, Random design)",
	method = function(model, draw) {
	require(glmnet)


	# Take out first row and column of Sigma (corresponding to latent variaable 
	# Z) that no longer needs to be there (again, assuming one block)
    if(model$nblocks==1 & model$sig_blocks ==1){
        Sigma <- model$Sigma[2:nrow(model$Sigma), 2:ncol(model$Sigma)]
    } else{
    	stop("!(model$nblocks==1 & model$sig_blocks ==1)")
    }
	# Set R to exactly reflect actual (known) clusters
	R <- (Sigma !=0)*1

	# Construct design matrix: replace cluster with average of features,
	# leave other features alone

	# # Identify clustered features: row i in R contains all the features
 #    # that are in a cluster with feature i, so come up with a list containing
 #    # all the clusters and then remove the repeats
 #    # print("R:")
 #    # print(R)
 #    clusters <- list()
 #    for(i in 1:nrow(R)){
 #        clusters[[i]] <- which(R[i, ] > 0)
 #    }

 #    clusters <- unique(clusters)
 #    # Only care about clusters with more than one element (only ones that need
 #    # to be treated differently)
 #    # keep track of whether there's more than one cluster or not
 #    multiple <- sum(lengths(clusters) > 1) > 1

 	fit <- protolasso(x=draw$X, y=draw$y, R)

	return(list(selected_sets=fit$selected_sets,
		non_cluster_feats=fit$non_cluster_feats, beta=fit$beta))	
	}
)

lassoMB_phat_ideal <- new_method("lassoMB_phat_ideal",
	"Ideal Lasso SS (M&B) p_hats",
	method = function(model, draw, cutoff, q, pfer, B) {
		# Set R to exactly reflect actual (known) clusters
		R <- (model$Sigma !=0)*1
		# Check if used a model where signal cluster variable was removed; if
		# so, remove first row and column of Sigma (and therefore R)
		# print(str(R))
		# print(str(model$x))
		# if(ncol(R) > ncol(model$x)){
		# 	R <- R[2:nrow(R), 2:ncol(R)]
		# }
	fit <- stabselGreg(x=model$x, y=draw, fitfun=glmnet.lasso,
		cutoff=cutoff, q=q, B=B, sampling_type="MB"
		, R=R
		)
	# return(list(phat=fit$phat[, ncol(fit$phat)], beta=NA))
	return(list(phat=fit$feat_sel_props, beta=NA, R=R, selected=fit$selected))
	},
	# Cutoff probability (proportion of time a feature must be selected by base
	# feature selection method to be included by stability selection) and Per 
	# family error rate (expected number of false selections) for stability
	# selection
	settings = list(cutoff = 0.52, q=11, B=100)
)

lassoSS_phat_cor <- new_method("lassoSS_phat_cor",
	"Corr Lasso SS (S&S) p_hats",
	method = function(model, draw, cutoff, q, pfer) {
		# Set R to be the absolute value of the correlation matrix
		R <- abs(cor(model$x))
		fit <- stabselGreg(x=model$x, y=draw, fitfun=glmnet.lasso,
			cutoff=cutoff, q=q, B=1000, sampling_type="SS"
			, R=R
			)
		# return(list(phat=fit$phat[, ncol(fit$phat)], beta=NA))
		return(list(phat=fit$feat_sel_props, beta=NA))
	},
	# Cutoff probability (proportion of time a feature must be selected by base
	# feature selection method to be included by stability selection) and Per 
	# family error rate (expected number of false selections) for stability
	# selection
	settings = list(cutoff = 0.7, q=11)
)

lassoSS_phat_cor_squared <- new_method("lassoSS_phat_cor_squared",
	"Corr Sq Lasso SS (S&S) p_hats",
	method = function(model, draw, cutoff, q, pfer, B) {
		# Set R to be the square of the correlation matrix
		R <- cor(model$x)^2
		fit <- stabselGreg(x=model$x, y=draw, fitfun=glmnet.lasso,
			cutoff=cutoff, q=q, B=B, sampling_type="SS"
			, R=R
			)
		# return(list(phat=fit$phat[, ncol(fit$phat)], beta=NA))
		return(list(phat=fit$feat_sel_props, beta=NA, B=B))
	},
	# Cutoff probability (proportion of time a feature must be selected by base
	# feature selection method to be included by stability selection) and Per 
	# family error rate (expected number of false selections) for stability
	# selection
	settings = list(cutoff = 0.7, q=11, B=100)
)

lassoMB_phat_cor_squared <- new_method("lassoMB_phat_cor_squared",
	"Corr Sq Lasso SS (M&B) p_hats",
	method = function(model, draw, cutoff, q, pfer, B) {
		# Set R to be the square of the correlation matrix
		R <- cor(model$x)^2
		fit <- stabselGreg(x=model$x, y=draw, fitfun=glmnet.lasso,
			cutoff=cutoff, q=q, B=B, sampling_type="MB"
			, R=R
			)
		# return(list(phat=fit$phat[, ncol(fit$phat)], beta=NA))
		return(list(phat=fit$feat_sel_props, beta=NA))
	},
	# Cutoff probability (proportion of time a feature must be selected by base
	# feature selection method to be included by stability selection) and Per 
	# family error rate (expected number of false selections) for stability
	# selection
	settings = list(cutoff = 0.7, q=11, B=100)
)

# #

# lasso_ss_pars_.9b <- matrix(c(0.9, 1500, 0.9, 1900),
# 	ncol=2, byrow=T)

# list_of_lasso_sss_5_.9b <- make_lasso_ss_columns(lasso_ss_pars_.9b)

# lasso_ss_pars_.7b <- matrix(c(0.7, 1000, 0.7, 1500, 0.7, 1900),
# 	ncol=2, byrow=T)

# list_of_lasso_sss_5_.7b <- make_lasso_ss_columns(lasso_ss_pars_.7b)

# lasso_ss_pars_.51b <- matrix(c(0.51, 250, 0.51, 500), ncol=2,
# 	byrow=T)

# list_of_lasso_sss_5_.51b <- make_lasso_ss_columns(lasso_ss_pars_.51b)

# ### Vanilla (cross-validated) lasso

# cv.lasso <- new_method("cv.lasso", "CV_Lasso",
# 	method = function(model, draw) {
# 	fit <-cv.glmnet(x=model$x, y=draw, family="gaussian", alpha=1, parallel=T)
# 	coefs <- coef(fit, s=fit$lambda.min)
# 	selected <- which(coefs[-1]!=0)
# 	return(list(selected=selected))
# 	}
# )

# ### Vanilla lasso (with q set as in stability selection)

# lasso <- new_method("lasso", "Lasso",
# 	method = function(model, draw, q) {
# 	fit <- suppressWarnings(glmnet(x=model$x, y=draw, family="gaussian",
# 		alpha=1, pmax=q))
# 	selected <- predict(fit, type="nonzero")
#     # Choose the features chosen at the particular lambda where q
#     # features were selected (i.e. the smallest or last lambda)
#     selected <- selected[[length(selected)]]
# 	return(list(selected=selected))
# 	}, 
# 	settings=list(q = 200)
# )

### Vanilla lasso

lasso <- new_method("lasso", "Lasso",
	# Currently with a large lambda.min.ratio to ensure smaller selected sets.
	method = function(model, draw) {
	fit <- glmnet(x=model$x, y=draw, family="gaussian", alpha=1,
		nlambda=2000, lambda.min.ratio=0.1)
	selected <- predict(fit, type="nonzero")
	return(list(selected=selected, beta=fit$beta))
	}
)

lasso_random <- new_method("lasso_random", "Lasso (Random design)",
	# Currently with a large lambda.min.ratio to ensure smaller selected sets.
	method = function(model, draw, model_size) {
	fit <- glmnet(x=draw$X, y=draw$y, family="gaussian", alpha=1,
		nlambda=2000, lambda.min.ratio=0.1)
	selected <- unique(predict(fit, type="nonzero"))

	cond <- is.null(selected[[1]])
	while(cond){
		selected <- selected[2:length(selected)]
		cond <- is.null(selected[[1]])
	}
	selected <- selected[lengths(selected) <= model_size]

	return(list(lasso_selected=selected, testX=draw$testX, testY=draw$testY,
		testMu=draw$testMu))
	},
	settings = list(model_size=11)
)

elastic_net <- new_method("elastic_net", "Elastic Net",
	# Currently with a large lambda.min.ratio to ensure smaller selected sets.
	method = function(model, draw, model_size) {
	fit <- glmnet(x=draw$X, y=draw$y, family="gaussian", alpha=0.5,
		nlambda=2000, lambda.min.ratio=0.1)
	selected <- unique(predict(fit, type="nonzero"))

	cond <- is.null(selected[[1]])
	while(cond){
		selected <- selected[2:length(selected)]
		cond <- is.null(selected[[1]])
	}
	selected <- selected[lengths(selected) <= model_size]

	return(list(lasso_selected=selected, testX=draw$testX, testY=draw$testY,
		testMu=draw$testMu))
	},
	settings = list(model_size=11)
)

### Lasso of size 2

lasso2 <- new_method("lasso2", "Lasso2",
	method = function(model, draw) {
		require(purrr)
		require(glmnet)
		# get correlation matrix
		cor_mat <- cor(cbind(model$x, draw))
		n_lambda <- ncol(model$x)*10
		while(n_lambda > 0){
			fit <- glmnet(x=model$x, y=draw, family="gaussian",
				nlambda=n_lambda, alpha=1)
			selected_sets <- unique(predict(fit, type="nonzero")) %>%
				discard(is.null)
			for(i in 1:length(selected_sets)){
				if(length(selected_sets[[i]]) == 2){
					# print("selected:")
					# print(selected_sets[[i]])
					# print("correlation matrix:")
					# print(cor_mat)
					return(list(selected=selected_sets[[i]], cor_mat=cor_mat))
				}
			}
		# if wasn't able to return a value, try again with bigger n_lambda
		n_lambda <- n_lambda*10
		}
		
	}
)

# Lasso fit on sample of size floor(n/2)

lasso.n2 <- new_method("lasso.n2", "Lasso_n2",
	method = function(model, draw) {
		n <- length(draw)
		samp_size <- floor(n/2)
		order <- sample(n)[1:samp_size]
		fit <- glmnet(x=model$x[order, ], y=draw[order], family="gaussian",
			alpha=1)
		selected <- predict(fit, type="nonzero")
		return(list(selected=selected, beta=fit$beta))
	}
)

# make_lasso <- function(qs) {
# 	# qs is an n-vector. The ith entry is the desired number of features
# 	# selected in each lasso fit (q) 
# 	ret <- list()
# 	for(i in 1:length(qs)){
# 		meth <- new_method(name = sprintf("lasso_%s", qs[i]),
#         label = sprintf("Lasso (q=%s)", qs[i]),
#         method = function(model, draw, q) {
# 			fit <- suppressWarnings(glmnet(x=model$x, y=draw, family="gaussian",
# 				alpha=1, pmax=q))
# 			selected <- predict(fit, type="nonzero")
# 			# Choose the features chosen at the particular lambda where q
# 		    # features were selected (i.e. the smallest or last lambda)
# 		    selected <- selected[[length(selected)]]
# 			return(list(selected=selected))
# 		}, 
# 		# Cutoff probability (proportion of time a feature must be selected by base
# 		# feature selection method to be included by stability selection) and Per 
# 		# family error rate (expected number of false selections) for stability
# 		# selection
# 		settings = list(q = qs[i])
# 		)
# 		ret <- c(ret, meth)
# 	}
#   	return(ret)
# }

# list_of_lassos <- make_lasso(qs)

# make_lasso <- function(q) {
# 	meth <- new_method(name = sprintf("lasso_%s", q),
#     label = sprintf("Lasso (q=%s)", q),
#     method = function(model, draw, q) {
# 		fit <- suppressWarnings(glmnet(x=model$x, y=draw, family="gaussian",
# 			alpha=1, pmax=q))
# 		selected <- predict(fit, type="nonzero")
# 		# Choose the features chosen at the particular lambda where q
# 	    # features were selected (i.e. the smallest or last lambda)
# 	    selected <- selected[[length(selected)]]
# 		return(list(selected=selected))
# 	}, 
# 	# Cutoff probability (proportion of time a feature must be selected by base
# 	# feature selection method to be included by stability selection) and Per 
# 	# family error rate (expected number of false selections) for stability
# 	# selection
# 	settings = list(q = q)
# 	)
# 	ret <- c(ret, meth)
#   	return(ret)
# }


### (Cross-validated) Elastic Net

# cv.elastic.net <- new_method("cv.elastic.net", "CV_EN",
# 	method = function(model, draw, alpha) {
# 	## Select best alpha
# 	# Sequence of alphas to use (concentrated near 0)
# 	alphas <- seq(0, 1, len = 11)^3
# 	search <- foreach(i = alphas, .combine = rbind) %dopar% {
# 		cv <- cv.glmnet(x=model$x, y=draw, family = "gaussian", alpha = i,
# 			parallel=T)
# 		data.frame(cvm = cv$cvm[cv$lambda == cv$lambda.min],
# 			lambda.min = cv$lambda.min, alpha = i)
# 	}
# 	# Identify lowest alpha and lambda combination
# 	cv3 <- search[search$cvm == min(search$cvm), ]

# 	# Fit using this combination
# 	fit <- glmnet(x=model$x, y=draw, family = "gaussian",
# 		lambda = cv3$lambda.min, alpha = cv3$alpha)
# 	coefs <- coef(fit)
# 	# fit <-cv.glmnet(x=model$x, y=draw, family="gaussian",
# 	# 	alpha=alpha)
# 	# fit <-cva.glmnet(x=model$x, y=draw, family="gaussian")
# 	# coefs <- coef(fit, s=fit$lambda.min)
# 	selected <- which(coefs[-1]!=0)
# 	# print("Elastic Net Selected features:")
# 	# print(selected)
# 	# print(typeof(selected))
# 	return(list(selected=selected))
# 	}
# 	# },
# ## alpha for elastic net (later choose using cross validation, maybe with this
# ## package: https://github.com/hong-revo/glmnetUtils)
# 	# settings = list(alpha = 0.5)
# )

# cv.elastic.net.alpha.0.5 <- new_method("cv.elastic.net.alpha.0.5", "CV_EN_0.5",
# 	method = function(model, draw) {
# 		cv <- cv.glmnet(x=model$x, y=draw, family = "gaussian", alpha = 0.5)
# 		# Fit using best lambda combination
# 		fit <- glmnet(x=model$x, y=draw, family = "gaussian", alpha = 0.5,
# 			lambda = cv$lambda.min)
# 		coefs <- coef(fit)
# 		selected <- which(coefs[-1]!=0)
# 		return(list(selected=selected))
# 	}
# )

# ### alpha=0.5 Elastic net (with q set as in stability selection)

# elastic.net.alpha.0.5 <- new_method("elastic.net.alpha.0.5", "EN_0.5",
# 	method = function(model, draw, q) {
# 	fit <- suppressWarnings(glmnet(x=model$x, y=draw, family="gaussian",
# 		alpha=0.5, pmax=q))
# 	selected <- predict(fit, type="nonzero")
#     # Choose the features chosen at the particular lambda where q
#     # features were selected (i.e. the smallest or last lambda)
#     selected <- selected[[length(selected)]]
# 	return(list(selected=selected))
# 	}, 
# 	settings=list(q = 200)
# )

### alpha=0.5 Elastic net

elastic.net.alpha.0.5 <- new_method("elastic.net.alpha.0.5", "EN_0.5",
	method = function(model, draw) {
	fit <- glmnet(x=model$x, y=draw, family="gaussian", alpha=0.5)
	selected <- predict(fit, type="nonzero")
    # Choose the features chosen at the particular lambda where q
    # features were selected (i.e. the smallest or last lambda)
	return(list(selected=selected, beta=fit$beta))
	}
)

# Subspace Stability Selection (proposed)
subspace_ss_max_cancor <- new_method(name = "subspace_ss_max_cancor",
	label = "Subspace_SS_Max_Cancor",
	method = function(model, draw, B, q) {
		require(ccaPP)
		require(Matrix)
		n <- nrow(model$x)
		p <- ncol(model$x)
		if(is.na(q)){
			q <- min(5, round(0.5*(model$sig_blocks + model$k_unblocked)))
			q <- max(2, q)
		}
		# Select 2*B iterations of features by complementary pair sampling
		selected_features <- complementaryPairSamplingIndicators(X.dat=model$x,
			y=draw, B=B, q=q)
		# Now selected_features contains 2*B columns of selected features.
		# For each set of selected features, calculate all of the the maximum
		# canonical correlations between the selected feature set and the
		# remaining selected feature sets. Take the median of these canonical
		# correlations as our metric for closeness to the center.

		# # (2*B - 1) by 2*B matrix; each of the 2*B columns will contain the
		# # (2*B - 1) maximum canonical correlations the ith feature space has
		# # with each of the (2*B - 1) other feature spaces. (I will skip over the
		# # canonical correlation each feature space has with itelf.)
		# cancors <- matrix(numeric(2*B*(2*B-1)), ncol=(2*B))
		# # 2*B-dimensional vector that will contain the median maximum canonical
		# # correlation each of the 2*B feature spaces has with the other 2*B - 1
		# # feature spaces.
		# median_cancors <- numeric(2*B)
		# for(i in 1:(2*B)){
		# 	# feature sets I will iterate over (all but feature set i)
		# 	j_feats <- (1:(2*B))[-i]
		# 	for(j in 1:(2*B-1)){
		# 		# canonical correlations
		# 		# cancor_output <- cancor(model$x %*% selected_features[, i],
		# 		# 	model$x %*% selected_features[, j_feats[j]])
		# 		# # count maximum canonical correlation between selected
		# 		# # feature sets i and j.
		# 		# cancors[j, i] <- max(cancor_output[[1]])
		# 		cancors[j, i] <- ccaPP::maxCorGrid(model$x[,
		# 		 	 	selected_features[, i]], model$x[,
		# 		 	 	selected_features[, j_feats[j]]], method="pearson")$cor
		# 	}
		# }
		# # Make medican_cancors contain the median of each column of canonical
		# # correlations.
		# median_cancors <- apply(cancors, 2, median)


		# 2*B by 2*B matrix; entry ij contains the canonical correlation between
		# feature spaces i and j.
		cancors <- matrix(numeric(2*B*2*B), ncol=(2*B))
		# 2*B-dimensional vector that will contain the median maximum canonical
		# correlation each of the 2*B feature spaces has with the other 2*B - 1
		# feature spaces.
		median_cancors <- numeric(2*B)
		# Calculate QR decompositions of 2*B matrices in advance
		qrs <- list()
		for(i in 1:(2*B)){
			x_i <- as.matrix(model$x[, as.logical(selected_features[, i])])
			qr_i <- qr(x_i)
			qrs[[i]] <- qr_i
			# if(length(dim(x_i)) == 1){
			# 	x_i <- as.matrix(x_i)
			# 	rank_i <- 1
			# } else{
				rank_i <- rankMatrix(x_i)
			# }
			# Check if rank makes sense
			if(rank_i != qr_i$rank){
				print(rankMatrix(x_i))
				print(qr_i$rank)
				stop("rankMatrix(x_i) != qr_i$rank")
			}
			# if(rank_i==1){
			# 	print("rank_i == 1")
			# }
		}
		for(i in 1:(2*B-1)){
			for(j in (i+1):(2*B)){
				# canonical correlations
				cancor_ij <- canCorBase(model$x,
					as.logical(selected_features[, i]),
					as.logical(selected_features[, j]), XA_QR=qrs[[i]],
					XB_QR=qrs[[j]], max=TRUE)
				cancors[i, j] <- cancor_ij
				cancors[j, i] <- cancor_ij
			}
		}
		# Make medican_cancors contain the median of each column of canonical
		# correlations (except canonical correlation of space with itself)
		for(i in 1:(2*B)){
			cors_i <- cancors[, i]
			# Remove canonical correlation of space with itself
			cors_i <- cors_i[-i]
			median_cancors[i] <- median(cors_i)
		}
		# # The feature set with the largest median canonical correlation
		# # is the one we want.
		# stability_selected_indices <- median_cancors == max(median_cancors)
		# print(min(which(stability_selected_indices)))
		# # If there is a tie between highest median canonical correlation,
		# # arbitratily choose the first feature space.
		# stability_selected_index <- min(which(stability_selected_indices))


		# Select the feature set with the largest median canonical correlation
		# (if there is a tie, pick the sparser set; if there is still a tie,
		# arbitrarily choose the first set)
		max <- which(median_cancors == max(median_cancors))
		if(length(max) > 1){
			# In event of tie, find sparsest selected set among ties
			sparsities <- integer(length(max))
			for(i in 1:length(max)){
				sparsities[i] <- sum(selected_features[, max[i]])
			}
			sparsest <- which(sparsities == min(sparsities))
			max <- max[sparsest]
			# If still tied, arbitrarily choose first selected set
			max <- max[1]
		}

		stability_selected  <- which(selected_features[, max] !=0)
		return(list(selected=stability_selected, beta=NA, q=q))
	},
	settings=list(B = 100, q = NA)
)

# Subspace Stability Selection (proposed)
subspace_ss_min_cancor <- new_method(name = "subspace_ss_min_cancor",
	label = "Subspace_SS_Min_Cancor",
	method = function(model, draw, B, q) {
		require(ccaPP)
		require(Matrix)
		n <- nrow(model$x)
		p <- ncol(model$x)
		if(is.na(q)){
			q <- min(5, round(0.5*(model$sig_blocks + model$k_unblocked)))
			q <- max(2, q)
		}
		# Select 2*B iterations of features by complementary pair sampling
		selected_features <- complementaryPairSamplingIndicators(X.dat=model$x,
			y=draw, B=B, q=q)
		# Now selected_features contains 2*B columns of selected features.
		# For each set of selected features, calculate all of the the maximum
		# canonical correlations between the selected feature set and the
		# remaining selected feature sets. Take the median of these canonical
		# correlations as our metric for closeness to the center.

		# # (2*B - 1) by 2*B matrix; each of the 2*B columns will contain the
		# # (2*B - 1) maximum canonical correlations the ith feature space has
		# # with each of the (2*B - 1) other feature spaces. (I will skip over the
		# # canonical correlation each feature space has with itelf.)
		# cancors <- matrix(numeric(2*B*(2*B-1)), ncol=(2*B))
		# # 2*B-dimensional vector that will contain the median maximum canonical
		# # correlation each of the 2*B feature spaces has with the other 2*B - 1
		# # feature spaces.
		# median_cancors <- numeric(2*B)
		# for(i in 1:(2*B)){
		# 	# feature sets I will iterate over (all but feature set i)
		# 	j_feats <- (1:(2*B))[-i]
		# 	for(j in 1:(2*B-1)){
		# 		# canonical correlations
		# 		# cancor_output <- cancor(model$x %*% selected_features[, i],
		# 		# 	model$x %*% selected_features[, j_feats[j]])
		# 		# # count maximum canonical correlation between selected
		# 		# # feature sets i and j.
		# 		# cancors[j, i] <- max(cancor_output[[1]])
		# 		cancors[j, i] <- ccaPP::maxCorGrid(model$x[,
		# 		 	 	selected_features[, i]], model$x[,
		# 		 	 	selected_features[, j_feats[j]]], method="pearson")$cor
		# 	}
		# }
		# # Make medican_cancors contain the median of each column of canonical
		# # correlations.
		# median_cancors <- apply(cancors, 2, median)


		# 2*B by 2*B matrix; entry ij contains the canonical correlation between
		# feature spaces i and j.
		cancors <- matrix(numeric(2*B*2*B), ncol=(2*B))
		# 2*B-dimensional vector that will contain the median maximum canonical
		# correlation each of the 2*B feature spaces has with the other 2*B - 1
		# feature spaces.
		median_cancors <- numeric(2*B)
		# Calculate QR decompositions of 2*B matrices in advance
		qrs <- list()
		for(i in 1:(2*B)){
			x_i <- as.matrix(model$x[, as.logical(selected_features[, i])])
			qr_i <- qr(x_i)
			qrs[[i]] <- qr_i
			# if(length(dim(x_i)) == 1){
			# 	x_i <- as.matrix(x_i)
			# 	rank_i <- 1
			# } else{
				rank_i <- rankMatrix(x_i)
			# }
			# Check if rank makes sense
			if(rank_i != qr_i$rank){
				print(rankMatrix(x_i))
				print(qr_i$rank)
				stop("rankMatrix(x_i) != qr_i$rank")
			}
			# if(rank_i==1){
			# 	print("rank_i == 1")
			# }
		}
		for(i in 1:(2*B-1)){
			for(j in (i+1):(2*B)){
				# canonical correlations
				cancor_ij <- canCorBase(model$x,
					as.logical(selected_features[, i]),
					as.logical(selected_features[, j]), XA_QR=qrs[[i]],
					XB_QR=qrs[[j]], max=FALSE)
				cancors[i, j] <- cancor_ij
				cancors[j, i] <- cancor_ij
			}
		}
		# Make medican_cancors contain the median of each column of canonical
		# correlations (except canonical correlation of space with itself)
		for(i in 1:(2*B)){
			cors_i <- cancors[, i]
			# Remove canonical correlation of space with itself
			cors_i <- cors_i[-i]
			median_cancors[i] <- median(cors_i)
		}
		# # The feature set with the largest median canonical correlation
		# # is the one we want.
		# stability_selected_indices <- median_cancors == max(median_cancors)
		# print(min(which(stability_selected_indices)))
		# # If there is a tie between highest median canonical correlation,
		# # arbitratily choose the first feature space.
		# stability_selected_index <- min(which(stability_selected_indices))


		# Select the feature set with the largest median canonical correlation
		# (if there is a tie, pick the sparser set; if there is still a tie,
		# arbitrarily choose the first set)
		max <- which(median_cancors == max(median_cancors))
		if(length(max) > 1){
			# In event of tie, find sparsest selected set among ties
			sparsities <- integer(length(max))
			for(i in 1:length(max)){
				sparsities[i] <- sum(selected_features[, max[i]])
			}
			sparsest <- which(sparsities == min(sparsities))
			max <- max[sparsest]
			# If still tied, arbitrarily choose first selected set
			max <- max[1]
		}

		stability_selected  <- which(selected_features[, max] !=0)
		return(list(selected=stability_selected, beta=NA, q=q))
	},
	settings=list(B = 100, q = NA)
)

# Canonical correlation
canCor <- function(X, selected_a, selected_b){
	# Inputs:
	#
	# X: data matrix
	# selected_a: first selected set
	# selected_b: second selected set
	#
	# Output: canoncial correlation
	require(ccaPP)
	# Space a
	X_a <- X[, selected_a]
	# Space b
	X_b <- X[, selected_b]
	# Calculate distance metric
	return(ccaPP::maxCorGrid(X_a, X_b, method="pearson")$cor)
}

# Canonical correlation
canCorBase <- function(X, selected_a, selected_b, XA_QR=NA, XB_QR=NA,
	SVD_AB=NA, max=TRUE){
	# Inputs:
	#
	# X: data matrix
	# selected_a: first selected set
	# selected_b: second selected set
	#
	# Output: canoncial correlation
	require(irlba)
	# Space a
	X_a <- X[, selected_a]
	# Space b
	X_b <- X[, selected_b]
	# Calculate distance metric
	# return(ccaPP::maxCorGrid(X_a, X_b, method="pearson")$cor)
	return(myCancor(X_a, X_b, max=max, XA_QR=XA_QR, XB_QR=XB_QR))	
}

myCancor <- function (x, y, max, xcenter=TRUE, ycenter=TRUE, XA_QR=NA,
	XB_QR=NA){
	require(irlba)
    mat_to_svd <- getCancorMat(x, y, xcenter, ycenter, XA_QR, XB_QR)
	return(getSV(mat_to_svd, max))
}

getCancorMat <- function(x, y, xcenter, ycenter, XA_QR, XB_QR){
	x <- as.matrix(x)
    y <- as.matrix(y)
    if ((nr <- nrow(x)) != nrow(y)) 
        stop("unequal number of rows in 'cancor'")
    ncx <- ncol(x)
    ncy <- ncol(y)
    if (!nr || !ncx || !ncy) 
        stop("dimension 0 in 'x' or 'y'")
    if (is.logical(xcenter)) {
        if (xcenter) {
            xcenter <- colMeans(x, )
            x <- x - rep(xcenter, rep.int(nr, ncx))
        }
        else xcenter <- rep.int(0, ncx)
    }
    else {
        xcenter <- rep_len(xcenter, ncx)
        x <- x - rep(xcenter, rep.int(nr, ncx))
    }
    if (is.logical(ycenter)) {
        if (ycenter) {
            ycenter <- colMeans(y)
            y <- y - rep(ycenter, rep.int(nr, ncy))
        }
        else ycenter <- rep.int(0, ncy)
    }
    else {
        ycenter <- rep_len(ycenter, ncy)
        y <- y - rep(ycenter, rep.int(nr, ncy))
    }
    # Calculate QR decompositions inverses
	if(length(XA_QR) == 1){
		if(is.na(XA_QR)){
			XA_QR <- qr(x)
		}
	}
	if(length(XB_QR) == 1){
		if(is.na(XB_QR)){
			XB_QR <- qr(y)
		}
	}
    qx <- XA_QR
    qy <- XB_QR
    dx <- qx$rank
    if (!dx) {
        stop("'x' has rank 0")
    }
    dy <- qy$rank
    if (!dy) {
        stop("'y' has rank 0")
    }
    # Matrix we need largest (or smallest) singular value of
    mat_to_svd <- as.matrix(qr.qty(qx, qr.qy(qy, diag(1, nr, dy)))[1L:dx, , 
        drop = FALSE])
    return(mat_to_svd)
}

getHashResult <- function(hash, key){
	res <- hash[[key]]
	if(is.null(res)){
		stop("Error: no response available for key.")
	}
	return(res)
}

# Subspace selection (Proposed 10/16/20)
subspace_sel_grass <- new_method(name = "subspace_sel_grass",
	label = "Subspace_Sel_Grass",
	method = function(model, draw, B, q) {
		# require(MASS)
		# require(hash)
		# require(digest)
		n <- nrow(model$x)
		p <- ncol(model$x)
		if(is.na(q)){
			q <- min(5, round(0.5*(model$sig_blocks + model$k_unblocked)))
			q <- max(2, q)
		}
		# Select 2*B iterations of features by complementary pair sampling
		selected_features <- complementaryPairSamplingIndicators(X.dat=model$x,
			y=draw, B=B, q=q)
		# Now selected_features contains 2*B columns of selected features.
		# For each pair of selected features, calculate the distance between
		# the pair according to the metric defined in the subDist function.

		# 2*B by 2*B matrix; entry (i, j) will contain the distance between
		# selected sets i and j.
		distances <- matrix(numeric(2*B*2*B), ncol=(2*B))

		# First calculate all 2*B Moore-Penrose inverses of the feature spaces
		# (will save time to calculate ahead of time). Will store in a hash
		# table with key generated from selected feature set in order to 
		# avoid re-calculating inverse of any matrix we've already calculated
		# MP_invs_hash <- hash()
		MP_invs <- list()
		for(i in 1:(2*B)){
			if(all(selected_features[, i] == 0)){
				print(i)
				print(dim(model$x[, selected_features[, i]]))
				print(selected_features[, i])
			}
			if(any(dim(model$x[, selected_features[, i]]) == 0)){
				print(i)
				print(dim(model$x[, selected_features[, i]]))
				print(selected_features[, i])
			}
			# key_i <- digest(selected_features[, i])
			# if(is.null(MP_invs_hash[[key_i]])){
			# 	MP_invs_hash[[key_i]] <- ginv(model$x[, selected_features[, i]])
			# }
			MP_invs[[i]] <- ginv(model$x[, selected_features[, i]])
		}
		# Create similar distances hash table
		# distances_hash <- hash()
		for(i in 1:(2*B - 1)){
			feats_i <- selected_features[, i]
			# key_i <- digest(feats_i)
			# ginv_i <- getHashResult(MP_invs_hash, key_i)
			for(j in (i+1):(2*B)){
				feats_j <- selected_features[, j]
				# key_ij <- digest(list(feats_i, feats_j))
				# if(is.null(distances_hash[[key_ij]])){
				# 	key_j <- digest(feats_j)
				# 	ginv_j <- getHashResult(MP_invs_hash, key_j)
				# 	distance_ij <- subDist(model$x, as.logical(feats_i),
				# 		as.logical(feats_j), ginv_i, ginv_j)
				# 	distances_hash[[key_ij]] <- distance_ij
				# 	key_ji <- digest(list(feats_j, feats_i))
				# 	distances_hash[[key_ji]] <- distance_ij
				# }else{
				# 	distance_ij <- distances_hash[[key_ij]]
				# }
				distance_ij <- subDist(model$x, as.logical(feats_i),
						as.logical(feats_j), MP_invs[[i]], MP_invs[[j]])
				distances[i, j] <- distance_ij
				distances[j, i] <- distance_ij
 			}
		}
		# Make sure distances is symmetric as expected
		if(!all.equal(distances, t(distances))){
			stop("!all.equal(distances, t(distances))")
		}
		# Select the feature set with the smallest row sum in distances (if 
		# there is a tie, pick the sparser set; if there is still a tie,
		# arbitrarily choose the first set)
		sums <- rowSums(distances)
		# min <- which(sums == min(sums))
		min <- which.min(sums)
		if(length(min) > 1){
			# In event of tie, find sparsest selected set among ties
			sparsities <- integer(length(min))
			for(i in 1:length(min)){
				sparsities[i] <- sum(selected_features[, min[i]])
			}
			sparsest <- which(sparsities == min(sparsities))
			min <- min[sparsest]
			# If still tied, arbitrarily choose first selected set
			min <- min[1]
		}
		stability_selected  <- which(selected_features[, min] !=0)
		return(list(selected=stability_selected, beta=NA, q=q))
	},
	settings=list(B = 100, q = NA)
)

# Distance metric
subDist <- function(X, selected_a, selected_b, XA_inv=NA, XB_inv=NA){
	# Inputs:
	#
	# X: data matrix
	# selected_a: first selected set
	# selected_b: second selected set
	# XA_inv (optional): Moore-Penrose inverse of X[, selected_a] (can calculate
	# in advance and supply to function to save computational resources)
	# XB_inv (optional): Moore-Penrose inverse of X[, selected_a] (can calculate
	# in advance and supply to function to save computational resources)
	#
	# Output: distance metric proposed at 10/16/20 meeting
	require(irlba)
	require(MASS)
	# Calculate Moore-Penrose inverses
	if(length(XA_inv) == 1){
		if(is.na(XA_inv)){
			XA_inv <- ginv(X[, selected_a])
		}
	}
	if(length(XA_inv) == 1){
		if(is.na(XB_inv)){
			XB_inv <- ginv(X[, selected_b])
		}
	}
	# Calculate distance metric: operator norm (largest singular value) of
	# difference of projection matrices

	# Matrix we need largest singular value of
    mat_to_svd <- as.matrix(X[, selected_a] %*% XA_inv - X[, selected_b] %*%
    	XB_inv)
	# return(getSV(mat_to_svd, max=TRUE))
	return(norm(mat_to_svd, "2"))
}

getSV <- function(mat, max){
	# Return largest singular value of mat
	require(irlba)
    if(ncol(mat) == 1){
    	# If this is a one dimensional matrix, then the singular value is
    	# just the 2-norm.
    	ret <- sqrt(as.numeric(t(mat) %*% mat))
    } else if(ncol(mat) == 1){
    	ret <- sqrt(as.numeric(mat %*% t(mat)))
    } else{
    	if(!max){
    		# If getting mininum SV, check if matrix is rank deficient; if so,
    		# return 0
    		rank <- rankMatrix(mat)
    		if(rank < min(dim(mat))){
    			return(0)
    		}
    	}
    	if(max){
    		z <- irlba(mat, nv=1, smallest=(!max)) 	
		    if(length(z$d) > 1){
		    	stop("length(z$d) > 1")
		    }
		    ret <- z$d
	    }else{
	    	svd <- svd(mat)
	    	ret <- min(svd$d)
	    }
    }
	return(ret)
}


# make_subspace_ss <- function(qs) {
# 	# qs is an n-vector. The ith entry is the desired number of features
# 	# selected in each lasso fit (q) 
# 	ret <- list()
# 	for(i in 1:length(qs)){
# 		meth <- new_method(name = sprintf("subspace_ss_max_cancor_%s", qs[i]),
#         label = sprintf("SS_SS (q=%s)", qs[i]),
#         method = function(model, draw, B, q) {
# 			n <- nrow(model$x)
# 		p <- ncol(model$x)
# 		# Per Samworth & Shah, take a sample (without replacement)
# 		# of size floor(n/2)
# 		samp_size <- floor(n/2)
# 		# p by 2B matrix; in each of 2B columns, the jth row will equal 1
# 		# if that feature was selected, and zero otherwise
# 		selected_features <- matrix(numeric(p*2*B), ncol=(2*B))
# 		# Complete B repetitions of the following:
# 		for(i in 1:B){
# 			# Per Samworth & Shah, take complementary pairwise samples
# 			order <- sample(n)
# 			sample1 <- order[1:samp_size]
# 			sample2 <- order[(samp_size+1):(samp_size+samp_size)]
# 			# Fit the lasso using glmnet, choosing range of lambdas such that a
# 			# maximum of q features are selected.
# 			# use suppressWarnings because otherwise it gives this useless
# 			# error about q working properly
# 			fit1 <- suppressWarnings(glmnet(x=model$x[sample1, ],
# 				y=draw[sample1], family="gaussian", pmax=q))
# 			fit2 <- suppressWarnings(glmnet(x=model$x[sample2, ],
# 				y=draw[sample2], family="gaussian", pmax=q))
# 			## which coefficients are non-zero?
# 		    selected1 <- predict(fit1, type="nonzero")
# 		    selected2 <- predict(fit2, type="nonzero")
# 		    # Choose the features chosen at the particular lambda where q
# 		    # features were selected (i.e. the smallest or last lambda)
# 		    selected1 <- selected1[[length(selected1)]]
# 		    selected2 <- selected2[[length(selected2)]]
# 		    # In the (2i - 1)th column of selected_features, make the jth entry
# 		    # equal 1 if that feature was selected in selected1, and leave it
# 		    # as zero otherwise. Likewise with the (2i)th column and selected2.
# 		    selected_features[selected1, 2*i-1] <- 1
# 		    selected_features[selected2, 2*i] <- 1
# 		}
# 		# Now selected_features contains 2*B columns of selected features.
# 		# For each set of selected features, calculate all of the the maximum
# 		# canonical correlations between the selected feature set and the
# 		# remaining selected feature sets. Take the median of these canonical
# 		# correlations as our metric for closeness to the center.

# 		# (2*B - 1) by 2*B matrix; each of the 2*B columns will contain the
# 		# (2*B - 1) maximum canonical correlations the ith feature space has
# 		# with each of the (2*B - 1) other feature spaces. (I will skip over the
# 		# canonical correlation each feature space has with itelf.)
# 		cancors <- matrix(numeric(2*B*(2*B-1)), ncol=(2*B))
# 		# 2*B-dimensional vector that will contain the median maximum canonical
# 		# correlation each of the 2*B feature spaces has with the other 2*B - 1
# 		# feature spaces.
# 		median_cancors <- numeric(2*B)
# 		for(i in 1:(2*B)){
# 			# feature sets I will iterate over (all but feature set i)
# 			j_feats <- (1:(2*B))[-i]
# 			for(j in 1:(2*B-1)){
# 				# canonical correlations
# 				cancor_output <- cancor(model$x %*% selected_features[, i],
# 					model$x %*% selected_features[, j_feats[j]])
# 				# count maximum canonical correlation between selected
# 				# feature sets i and j.
# 				cancors[j, i] <- max(cancor_output[[1]])
# 			}
# 		}
# 		# Make medican_cancors contain the median of each column of canonical
# 		# correlations.
# 		median_cancors <- apply(cancors, 2, median)
# 		# The feature set with the largest median canonical correlation
# 		# is the one we want.
# 		stability_selected_indices <- median_cancors == max(median_cancors)
# 		print(min(which(stability_selected_indices)))
# 		# If there is a tie between highest median canonical correlation,
# 		# arbitratily choose the first feature space.
# 		stability_selected_index <- min(which(stability_selected_indices))
# 		stability_selected  <- which(selected_features[,
# 			stability_selected_index] !=0)
# 		return(list(selected=stability_selected, beta=NA))
# 		}, 
# 		# Cutoff probability (proportion of time a feature must be selected by base
# 		# feature selection method to be included by stability selection) and Per 
# 		# family error rate (expected number of false selections) for stability
# 		# selection
# 		settings = list(B = 50, q = qs[i])
# 		)
# 		ret <- c(ret, meth)
# 	}
#   	return(ret)
# }

lassoFS <- function(vec, n, p, q){
	# vec contains the subsample design matrix as well a the subsample response.
	# Fit a lasso model to the data and detect the nonzero coefficient
	# (selected) features. Return a vector of length p containing a 1
	# in the indices of selected features and a 0 elsewhere.

	# Unroll vec
	x <- matrix(vec[1:(n*p)], ncol=p)
	y <- vec[(n*p + 1):length(vec)]
	# Fit the lasso using glmnet, choosing range of lambdas such that a
	# maximum of q features are selected.
	# use suppressWarnings because otherwise it gives this useless
	# error about q working properly
	fit <- suppressWarnings(glmnet(x=x, y=y, family="gaussian", pmax=q))
	## which coefficients are non-zero?
	selected <- predict(fit, type="nonzero")
	ret <- numeric(p)
	ret[selected] <- 1
	return(ret)
}
	

make_subspace_ss <- function(qs) {
	require(ccaPP)
	# qs is an n-vector. The ith entry is the desired number of features
	# selected in each lasso fit (q) 
	ret <- list()
	for(i in 1:length(qs)){
		meth <- new_method(name = sprintf("subspace_ss_max_cancor_%s", qs[i]),
        label = sprintf("SS_SS (q=%s)", qs[i]),
        method = function(model, draw, B, q) {
			n <- nrow(model$x)
		p <- ncol(model$x)
		# Select 2*B iterations of features by complementary pair sampling
		selected_features <- complementaryPairSamplingIndicators(X.dat=model$x,
			y=draw, B=B, q=q)
		# Now selected_features contains 2*B columns of selected features.
		# For each set of selected features, calculate all of the the maximum
		# canonical correlations between the selected feature set and the
		# remaining selected feature sets. Take the median of these canonical
		# correlations as our metric for closeness to the center.

		# (2*B - 1) by 2*B matrix; each of the 2*B columns will contain the
		# (2*B - 1) maximum canonical correlations the ith feature space has
		# with each of the (2*B - 1) other feature spaces. (I will skip over the
		# canonical correlation each feature space has with itelf.)
		cancors <- matrix(numeric(2*B*(2*B-1)), ncol=(2*B))
		# 2*B-dimensional vector that will contain the median maximum canonical
		# correlation each of the 2*B feature spaces has with the other 2*B - 1
		# feature spaces.
		median_cancors <- numeric(2*B)
		for(i in 1:(2*B)){
			# feature sets I will iterate over (all but feature set i)
			j_feats <- (1:(2*B))[-i]
			for(j in 1:(2*B-1)){
				# canonical correlations
				# cancor_output <- cancor(model$x %*% selected_features[, i],
				# 	model$x %*% selected_features[, j_feats[j]])
				# # count maximum canonical correlation between selected
				# # feature sets i and j.
				# cancors[j, i] <- max(cancor_output[[1]])
				cancors[j, i] <- ccaPP::maxCorGrid(model$x[,
				 	selected_features[, i]], model$x[,
				 	selected_features[, j_feats[j]]], method="pearson")$cor

			}
		}
		# Make medican_cancors contain the median of each column of canonical
		# correlations.
		median_cancors <- apply(cancors, 2, median)
		# The feature set with the largest median canonical correlation
		# is the one we want.
		stability_selected_indices <- median_cancors == max(median_cancors)
		print(min(which(stability_selected_indices)))
		# If there is a tie between highest median canonical correlation,
		# arbitratily choose the first feature space.
		stability_selected_index <- min(which(stability_selected_indices))
		stability_selected  <- which(selected_features[,
			stability_selected_index] !=0)
		return(list(selected=stability_selected, beta=NA))
		}, 
		# Cutoff probability (proportion of time a feature must be selected by base
		# feature selection method to be included by stability selection) and Per 
		# family error rate (expected number of false selections) for stability
		# selection
		settings = list(B = 50, q = qs[i])
		)
		ret <- c(ret, meth)
	}
  	return(ret)
}

list_of_subspace_ss <- make_subspace_ss(qs)

lassoFSVector <- function(vec, n, p, q){
	# Like lassoFS, except return coefficients.
	# vec contains the subsample design matrix as well a the subsample response.
	# Fit a lasso model to the data and detect the nonzero coefficient
	# (selected) features. Return a vector of length p containing a 1
	# in the indices of selected features and a 0 elsewhere.

	# Unroll vec
	x <- matrix(vec[1:(n*p)], ncol=p)
	y <- vec[(n*p + 1):length(vec)]
	# Fit the lasso using glmnet, choosing range of lambdas such that a
	# maximum of q features are selected.
	# use suppressWarnings because otherwise it gives this useless
	# error about q working properly
	fit <- suppressWarnings(glmnet(x=x, y=y, family="gaussian", pmax=q))
	return(fit$beta[, ncol(fit$beta)])
}

# Vector Distance Stability Selection (proposed)
vector.d.ss <- new_method(name = "vector.d.ss",
	label = "Vector_D_SS",
	method = function(model, draw, B, q) {
		n <- nrow(model$x)
		p <- ncol(model$x)
		# Select 2*B iterations of features by complementary pair sampling
		feature_coefficients <- complementaryPairSamplingCoefficients(X.dat=
			model$x, y=draw, B=B, q=q)
		# Now feature_coefficients contains 2*B columns of coefficients.
		# For each set of coefficients beta_i, calculate all of the the 
		# Pearson correlations between X %*% beta_i and X %*% beta_j, where
		# beta_j iterates across the other (2*B - 1) lasso fits. Take the median
		# of these correlations as our metric for closeness to the center.

		# (2*B - 1) by 2*B matrix; each of the 2*B columns will contain the
		# Euclidean distance beween the ith vector X %*% beta_i and each of
		# the (2*B - 1) other vectors X %*% beta_j. (I will skip over the
		# distance between each vector and itelf.)
		distances <- matrix(numeric(2*B*(2*B-1)), ncol=(2*B))
		# 2*B-dimensional vector that will contain the median distance each
		# of the 2*B vectors has with the other 2*(B - 1) vectors.
		median_distances <- numeric(2*B)
		for(i in 1:(2*B)){
			# feature sets I will iterate over (all but feature set i)
			j_feats <- (1:(2*B))[-i]
			for(j in 1:(2*B-1)){
				# distances
				distances[j, i] <- as.numeric(dist(rbind(model$x %*% 
					feature_coefficients[, i], model$x %*% 
					feature_coefficients[, j_feats[j]])))
			}
		}
		# Make medican_distances contain the median of each column of correlations.
		median_distances <- apply(distances, 2, median)
		# The feature set with the smallest median distance
		# is the one we want.
		stability_selected_indices <- median_distances == min(median_distances)
		print(min(which(stability_selected_indices)))
		# If there is a tie between lowest median distance,
		# arbitratily choose the first feature space.
		stability_selected_index <- min(which(stability_selected_indices))
		stability_selected  <- which(feature_coefficients[,
			stability_selected_index] !=0)
		# Even though I have convenient access to the selected beta in
		# feature_coefficients[, stability_selected_index], it will be easier
		# for me in the evaluation stage if all of the stability selection
		# methods have beta = NA.
		return(list(selected=stability_selected, beta=NA))
	},
	settings=list(B = 50, q = 200)
)

# Create a matrix with 2B columns. Each column contains floor(n/2) X values
# followed by corresponding y values. Run an apply loop across these 2B columns
# using a function that recognizes which are the xs and which are the ys,
# then applies glmnet and returns results.

# vector_ss_function <- function(X.dat, y.dat){
# 	# Per Samworth & Shah, take complementary pairwise samples
# 		order <- sample(n)
# 		sample1 <- order[1:samp_size]
# 		sample2 <- order[(samp_size+1):(samp_size+samp_size)]
# 		# use suppressWarnings because otherwise it gives this useless
# 		# error about q working properly
# 		fit1 <- suppressWarnings(glmnet(x=X.dat[sample1, ],
# 			y=y.dat[sample1], family="gaussian", pmax=q))
# 		fit2 <- suppressWarnings(glmnet(x=X.dat[sample2, ],
# 			y=y.dat[sample2], family="gaussian", pmax=q))
# 	    # Make the (2*i-1)th column of feature_coefficients equal the
# 	    # coefficients from lasso fit 1 (choosing the coefficients when q
# 	    # features are selected), and the (2*i)th column has coefficients
# 	    # from lasso fit 2.
# 	    return(c(fit1$beta[, ncol(fit1$beta)], fit2$beta[, ncol(fit2$beta)]))
# }

make_vector_ss <- function(qs) {
	# qs is an n-vector. The ith entry is the desired number of features
	# selected in each lasso fit (q) 
	ret <- list()
	for(i in 1:length(qs)){
		meth <- new_method(name = sprintf("vector.ss_%s", qs[i]),
        label = sprintf("V_SS (q=%s)", qs[i]),
        method = function(model, draw, B, q) {
        n <- nrow(model$x)
		p <- ncol(model$x)
		# Select 2*B iterations of features by complementary pair sampling
		feature_coefficients <- complementaryPairSamplingCoefficients(X.dat=
			model$x, y=draw, B=B, q=q)
		# Now feature_coefficients contains 2*B columns of coefficients.
		# For each set of coefficients beta_i, calculate all of the the 
		# Pearson correlations between X %*% beta_i and X %*% beta_j, where
		# beta_j iterates across the other (2*B - 1) lasso fits. Take the median
		# of these correlations as our metric for closeness to the center.

		# (2*B - 1) by 2*B matrix; each of the 2*B columns will contain the
		# (2*B - 1) correlations the ith vector X %*% beta_i has with each of
		# the (2*B - 1) other vectors X %*% beta_j. (I will skip over the
		# correlation each vector has with itelf.)
		cors <- matrix(numeric(2*B*(2*B-1)), ncol=(2*B))
		# 2*B-dimensional vector that will contain the median correlation each
		# of the 2*B vectors has with the other 2*(B - 1) vectors.
		median_cors <- numeric(2*B)
		for(i in 1:(2*B)){
			# feature sets I will iterate over (all but feature set i)
			j_feats <- (1:(2*B))[-i]
			for(j in 1:(2*B-1)){
				# canonical correlations
				cors[j, i] <- cor(model$x %*% feature_coefficients[, i],
					model$x %*% feature_coefficients[, j_feats[j]])
			}
		}
		# Make medican_cors contain the median of each column of correlations.
		median_cors <- apply(cors, 2, median)
		# The feature set with the largest median canonical correlation
		# is the one we want.
		stability_selected_indices <- median_cors == max(median_cors)
		print(min(which(stability_selected_indices)))
		# If there is a tie between highest median canonical correlation,
		# arbitratily choose the first feature space.
		stability_selected_index <- min(which(stability_selected_indices))
		stability_selected  <- which(feature_coefficients[,
			stability_selected_index] !=0)
		# Even though I have convenient access to the selected beta in
		# feature_coefficients[, stability_selected_index], it will be easier
		# for me in the evaluation stage if all of the stability selection
		# methods have beta = NA.
		return(list(selected=stability_selected, beta=NA))
		}, 
		# Cutoff probability (proportion of time a feature must be selected by base
		# feature selection method to be included by stability selection) and Per 
		# family error rate (expected number of false selections) for stability
		# selection
		settings = list(B = 50, q = qs[i])
		)
		ret <- c(ret, meth)
	}
  	return(ret)
}

list_of_vector_ss <- make_vector_ss(qs)

make_vector_D_ss <- function(qs) {
	# qs is an n-vector. The ith entry is the desired number of features
	# selected in each lasso fit (q) 
	ret <- list()
	for(i in 1:length(qs)){
		meth <- new_method(name = sprintf("vector.d.ss_%s", qs[i]),
        label = sprintf("V_D_SS (q=%s)", qs[i]),
        method = function(model, draw, B, q) {
        n <- nrow(model$x)
		p <- ncol(model$x)
		# Select 2*B iterations of features by complementary pair sampling
		feature_coefficients <- complementaryPairSamplingCoefficients(X.dat=
			model$x, y=draw, B=B, q=q)
		# Now feature_coefficients contains 2*B columns of coefficients.
		# For each set of coefficients beta_i, calculate all of the the 
		# Pearson correlations between X %*% beta_i and X %*% beta_j, where
		# beta_j iterates across the other (2*B - 1) lasso fits. Take the median
		# of these correlations as our metric for closeness to the center.

		# (2*B - 1) by 2*B matrix; each of the 2*B columns will contain the
		# Euclidean distance beween the ith vector X %*% beta_i and each of
		# the (2*B - 1) other vectors X %*% beta_j. (I will skip over the
		# distance between each vector and itelf.)
		distances <- matrix(numeric(2*B*(2*B-1)), ncol=(2*B))
		# 2*B-dimensional vector that will contain the median distance each
		# of the 2*B vectors has with the other 2*(B - 1) vectors.
		median_distances <- numeric(2*B)
		for(i in 1:(2*B)){
			# feature sets I will iterate over (all but feature set i)
			j_feats <- (1:(2*B))[-i]
			for(j in 1:(2*B-1)){
				# vectors to take distances of, arranged in a matrix
				vecs <- rbind(t(model$x %*% 
					feature_coefficients[, i]), t(model$x %*% 
					feature_coefficients[, j_feats[j]]))
				# print(dim(vecs))
				# distances
				distances[j, i] <- as.numeric(dist(vecs))
			}
		}
		# Make medican_distances contain the median of each column of correlations.
		median_distances <- apply(distances, 2, median)
		# The feature set with the smallest median distance
		# is the one we want.
		stability_selected_indices <- median_distances == min(median_distances)
		print(min(which(stability_selected_indices)))
		# If there is a tie between lowest median distance,
		# arbitratily choose the first feature space.
		stability_selected_index <- min(which(stability_selected_indices))
		stability_selected  <- which(feature_coefficients[,
			stability_selected_index] !=0)
		# Even though I have convenient access to the selected beta in
		# feature_coefficients[, stability_selected_index], it will be easier
		# for me in the evaluation stage if all of the stability selection
		# methods have beta = NA.
		return(list(selected=stability_selected, beta=NA))
		}, 
		# Cutoff probability (proportion of time a feature must be selected by base
		# feature selection method to be included by stability selection) and Per 
		# family error rate (expected number of false selections) for stability
		# selection
		settings = list(B = 50, q = qs[i])
		)
		ret <- c(ret, meth)
	}
  	return(ret)
}

list_of_vector_D_ss <- make_vector_D_ss(qs)

# Subspace stability selection with bootstrapped samples of size n

make_bs_subspace_ss <- function(qs) {
	require(ccaPP)
	# qs is an n-vector. The ith entry is the desired number of features
	# selected in each lasso fit (q) 
	ret <- list()
	for(i in 1:length(qs)){
		meth <- new_method(name = sprintf("subspace.bs_ss_%s", qs[i]),
        label = sprintf("SS_BS_SS (q=%s)", qs[i]),
        method = function(model, draw, B, q) {
		n <- nrow(model$x)
		p <- ncol(model$x)
		# Select B sets of features by bootstrap sampling of size n
		selected_features <- bootstrapSamplingIndicators(X.dat=model$x, y=draw,
			B=B, q=q)
		# Now selected_features contains B columns of selected features.
		# For each set of selected features, calculate all of the the maximum
		# canonical correlations between the selected feature set and the
		# remaining selected feature sets. Take the median of these canonical
		# correlations as our metric for closeness to the center.

		# (B - 1) by B matrix; each of the B columns will contain the
		# (B - 1) maximum canonical correlations the ith feature space has
		# with each of the (B - 1) other feature spaces. (I will skip over the
		# canonical correlation each feature space has with itelf.)
		cancors <- matrix(numeric(B*(B-1)), ncol=(B))
		# B-dimensional vector that will contain the median maximum canonical
		# correlation each of the B feature spaces has with the other B - 1
		# feature spaces.
		median_cancors <- numeric(B)
		for(i in 1:B){
			# feature sets I will iterate over (all but feature set i)
			j_feats <- (1:B)[-i]
			for(j in 1:(B-1)){
				if((!(qr(model$x %*% selected_features[, i]))$rank) |
					!(qr(model$x %*% selected_features[, j_feats[j]])$rank)) {
					cancors[j, i] <- 0
					print(paste("rank error on i = ", i, "j = ", j))
				} else{
					# canonical correlations
					# cancor_output <- cancor(model$x %*% selected_features[, i],
					# 	model$x %*% selected_features[, j_feats[j]])
					# # count maximum canonical correlation between selected
					# # feature sets i and j.
					# cancors[j, i] <- max(cancor_output[[1]])
					cancors[j, i] <- ccaPP::maxCorGrid(model$x[,
				 	 	selected_features[, i]], model$x[,
				 	 	selected_features[, j_feats[j]]], method="pearson")$cor
				}
				
			}
		}
		# Make medican_cancors contain the median of each column of canonical
		# correlations.
		median_cancors <- apply(cancors, 2, median)
		# The feature set with the largest median canonical correlation
		# is the one we want.
		stability_selected_indices <- median_cancors == max(median_cancors)
		print(min(which(stability_selected_indices)))
		# If there is a tie between highest median canonical correlation,
		# arbitratily choose the first feature space.
		stability_selected_index <- min(which(stability_selected_indices))
		stability_selected  <- which(selected_features[,
			stability_selected_index] !=0)
		return(list(selected=stability_selected, beta=NA))
		}, 
		# Cutoff probability (proportion of time a feature must be selected by base
		# feature selection method to be included by stability selection) and Per 
		# family error rate (expected number of false selections) for stability
		# selection
		settings = list(B = 50, q = qs[i])
		)
		ret <- c(ret, meth)
	}
  	return(ret)
}

list_of_bs_subspace_ss <- make_bs_subspace_ss(qs)

# Vector stability selection with bootstrapped samples of size n

make_bs_vector_ss <- function(qs) {
	# qs is an n-vector. The ith entry is the desired number of features
	# selected in each lasso fit (q) 
	ret <- list()
	for(i in 1:length(qs)){
		meth <- new_method(name = sprintf("vector.bs_ss_%s", qs[i]),
        label = sprintf("V_BS_SS (q=%s)", qs[i]),
        method = function(model, draw, B, q) {
        n <- nrow(model$x)
		p <- ncol(model$x)
		# Select B sets of features by bootstrap sampling of size n
		feature_coefficients <- bootstrapSamplingCoefficients(X.dat=model$x,
			y=draw,B=B, q=q)
		# Now feature_coefficients contains B columns of coefficients.
		# For each set of coefficients beta_i, calculate all of the the 
		# Pearson correlations between X %*% beta_i and X %*% beta_j, where
		# beta_j iterates across the other (B - 1) lasso fits. Take the median
		# of these correlations as our metric for closeness to the center.

		# (B - 1) by B matrix; each of the B columns will contain the
		# (B - 1) correlations the ith vector X %*% beta_i has with each of
		# the (B - 1) other vectors X %*% beta_j. (I will skip over the
		# correlation each vector has with itelf.)
		cors <- matrix(numeric(B*(B-1)), ncol=(B))
		# B-dimensional vector that will contain the median correlation each
		# of the B vectors has with the other (B - 1) vectors.
		median_cors <- numeric(B)
		for(i in 1:(B)){
			# feature sets I will iterate over (all but feature set i)
			j_feats <- (1:(B))[-i]
			for(j in 1:(B-1)){
				# canonical correlations
				cors[j, i] <- cor(model$x %*% feature_coefficients[, i],
					model$x %*% feature_coefficients[, j_feats[j]])
			}
		}
		# Make medican_cors contain the median of each column of correlations.
		median_cors <- apply(cors, 2, median)
		# The feature set with the largest median canonical correlation
		# is the one we want.
		stability_selected_indices <- median_cors == max(median_cors)
		print(min(which(stability_selected_indices)))
		# If there is a tie between highest median canonical correlation,
		# arbitratily choose the first feature space.
		stability_selected_index <- min(which(stability_selected_indices))
		stability_selected  <- which(feature_coefficients[,
			stability_selected_index] !=0)
		# Even though I have convenient access to the selected beta in
		# feature_coefficients[, stability_selected_index], it will be easier
		# for me in the evaluation stage if all of the stability selection
		# methods have beta = NA.
		return(list(selected=stability_selected, beta=NA))
		}, 
		# Cutoff probability (proportion of time a feature must be selected by base
		# feature selection method to be included by stability selection) and Per 
		# family error rate (expected number of false selections) for stability
		# selection
		settings = list(B = 50, q = qs[i])
		)
		ret <- c(ret, meth)
	}
  	return(ret)
}

list_of_bs_vector_ss <- make_bs_vector_ss(qs)

# Barber stability selection with bootstrapped samples of size n

# Auxilary functions:

getSparsePCs <- function(concat_mat, npc_max, len, colnames, p){
	# Selects sparsity tuning parameter and returns first npc_max sparse
	# principal components of concat_mat.

	# Inputs: concat_mat: a matrix
	# npc_max: the number of sparse principal components that should
	# be calculated.
	# len: the number of sparsity tuning parameters that should be 
	# considered in the cross-validation process to determine the best
	# sparsity tuning parameter.
	# colnames: a vector containing the column names of concat_mat, which are
	# simply numbers corresponding to the numbers of the features in the orignal
	# data matrix. This will allow me to correspond loadings in the sparse 
	# principal components to the features in the original design matrix.
	# p: the number of features in the original design matrix.

	# Output: spc.out.v, the sparse principal components, with dimension
	# p x npc_max. Each column of v contains a nonzero entry in row i if feature
	# i (in the original data matrix) has a nonzero loading in that 
	# principal component of matrix concat_mat).


	# Select sparsity tuning parameter using SPC.cv
	# print("finding sparsity tuning parameter...")
	bestsumabsv <- SPC.cv(concat_mat, trace=FALSE, sumabsv=seq(1.2,
		floor(sqrt(ncol(concat_mat))), len=len))$bestsumabsv
	# print("found!")
	# print(Sys.time() - t1)
	# Take the first npc_max principal components of this matrix.
	# print("finding principal components...")
	spc.out <- SPC(x=concat_mat, sumabsv=bestsumabsv, K=npc_max, trace=FALSE,
		compute.pve=FALSE)$v
	# Correspond these features to the original matrix.
	spc.out.v <- matrix(0, p, npc_max)
	for(j in 1:ncol(spc.out.v)){
		# Identify all the features of the original design matrix with nonzero
		# loadings in the jth column of spc.out.v
		# print("features:")
		# print(spc.out[, j] != 0)
		feats <- colnames[spc.out[, j] != 0]
		# print("also:")
		# print(feats)
		nonzero.loadings.j <- sort(unique(as.numeric(feats)))
		# Put ones in these entries of column j of spc.out.v
		spc.out.v[nonzero.loadings.j, j] <- 1
	}

	# print("done!")
	# print(Sys.time() - t1)
	return(spc.out.v)
}

getSelectedFeaturePathPC <- function(spc.out.v, npc_max, p){
	# Inputs: 
	# npc_max, the number of sparse principal components in
	# concat_mat
	# spc.out.v, a p x npc_max matrix containing the nonzero loadings
	# of the sparse principal components of concat_mat. Each column of 
	# spc.out.v is a factor (a sparse principal component).
	# p: number of features in the original design matrix
	#
	# Output: selected, a list of selected features of length npc_max.
	# Each selected feature set will be a vector, one for each number of 
	# principal components.

	selected <- list()
	# For each number of principal components, do the following:
	for(j in 1:npc_max){
		# Identify which features are nonzero in any one of the first j
		# principal components:
		# Identify all of the nonzero components in any of the columns,
		# one at a time.
		if(j==1){
			selected[[j]] <- (1:p)[spc.out.v[, 1]!= 0]
		}
		else{
			feat_sums <- rowSums(spc.out.v[, 1:j])
			selected[[j]] <- (1:p)[feat_sums != 0]
		}

		# # Add these selections to the list of selections
		# selected <- c(selected, selected_j)
	}
	return(selected)
}

# # pca proportions--proportions of q to use for npc
# pcaprops <- c(0.1, 0.2, 0.3, 0.4, 0.5)

# qs to use for barber subspace ss
qs_barber <- c(25, 50, 75, 100, 150, 200, 250, 300, 350, 400)

# Source code for PMA package: https://github.com/cran/PMA/blob/master/R/PMD.R

make_barber_subspace_ss3 <- function(qs_barber) {
	# qs_barber is an n-vector. The ith entry is the desired number of features
	# selected in each lasso fit (q) 
	ret <- list()
	for(i in 1:length(qs_barber)){
		meth <- new_method(name = sprintf("subspace.rb_ss3_%s", qs_barber[i]),
        label = sprintf("SS_RB_SS3 (q=%s)", qs_barber[i]),
        method = function(model, draw, B, q, npc_max) {
        t1 <- Sys.time()
		n <- nrow(model$x)
		p <- ncol(model$x)
		# Select B sets of features by bootstrap sampling of size n
		selected_features <- bootstrapSamplingIndicators(X.dat=model$x, y=draw,
			B=B, q=q)
		# print("selected_features:")
		# print(selected_features)
		# Now selected_features contains B columns of selected features.
		# Label the columns of model$x so that I will be able to identify
		# which features were selected after PCA.
		x.dat <- as.matrix(model$x)
		colnames(x.dat) <- as.character(1:p)
		# Form the concatenated n by sum(size(A_n)) matrix.
		concat_mat <- x.dat[, as.logical(selected_features[, 1])]
		for(i in 2:B){
			concat_mat <- cbind(concat_mat, x.dat[, selected_features[, i]])
		}
		# print("concat_mat:")
		# print(concat_mat[, 1:4])
		# print(dim(concat_mat))
		# Retrieve the first npc_max sparse principal componenets of concat_mat.
		spc.out.v <- getSparsePCs(concat_mat=concat_mat, npc_max=npc_max,
			len=20, colnames=colnames(concat_mat), p=p)
		# print("spc.out.v:")
		# print(spc.out.v)
		# Now spc.out.v is a p x npc_max matrix containing the nonzero loadings
		# of the sparse principal components of concat_mat in each column.

		# Next we will get a list of selected features, one at each number of 
		# principal components.
		selected <- getSelectedFeaturePathPC(spc.out.v=spc.out.v,
			npc_max=npc_max, p=p)
		# print("selected:")
		# print(selected)

		# Lastly, to make size eval work, create a matrix of betas which are
		# just 1 if the coefficient is nonzero and zero otherwise.
		beta <- matrix(0, p, npc_max)
		for(i in 1:npc_max){
			beta[selected[[i]], i] <- 1
		}
		# Return the selected features.
		# print("Complete")
		# print(Sys.time() - t1)
		return(list(selected=selected, beta=beta))
		}, 
		settings = list(B = 20, q = qs_barber[i], npc_max = 15)
		)
		ret <- c(ret, meth)
	}
  	return(ret)
}

list_of_barber_subspace_ss3 <- make_barber_subspace_ss3(c(10, 25, 50))

list_of_barber_subspace_ss3b <- make_barber_subspace_ss3(c(75, 100))

list_of_barber_subspace_ss3c <- make_barber_subspace_ss3(c(150, 200))

list_of_barber_subspace_ss3d <- make_barber_subspace_ss3(c(250, 300))

list_of_barber_subspace_ss3e <- make_barber_subspace_ss3(c(350, 400))

# # Function used for apply loop below.
# qrutil <- function(feats, X.dat){
# 	return(!(qr(X.dat[, feats==1])$rank))
# }

# q Max Cancor stability selection
#
# Fit B full lambda fits of lasso, using a sequence of nlambda lambdas. For each
# lambda, find the fit with the maximum canonical correlation with ALL other
# fits. Return that fit (for that lambda). Return a plot showing the tradeoff
# of screening criterion as the size of the selected set increases.


# Auxillary Functions

removeZeroLambdas <- function(selected_features, selected_features_lambdas,
	nlambda, B){
	# Output: 
	#
	# selected_features: a set of selected_features without any columns
	# corresponding to a lambda which included at least one selected feature set
	# with 0 features selected.
	#
	# selected_features_lambdas: a ncol(selected_features)-vector that
	# keeps track of the lambda corresponding to each selected feature set in
	# selected_features. The ith entry is the lambda corresponding to the
	# selected feature set in the ith column of selected_features. (This also
	# has the relevant lambdas removed, so that its length still matches the
	# number of columns of selected_features.)

	# Inputs: 
	#
	# selected_features: p x sum(n_feature_sets) matrix of selected features;
	# each column contains a selected features set, with a 1 in row i if feature
	# i was selected on that resample and lambda fit.
	#
	# selected_features_lambdas: a ncol(selected_features)-vector that
	# keeps track of the lambda corresponding to each selected feature set in
	# selected_features. The ith entry is the lambda corresponding to the
	# selected feature set in the ith column of selected_features.
	#
	# nlambda: number of lambdas used on lasso fits
	#
	# B: number of resamples used

	# Identify lambdas where an empty feature set was selected at least once:
	# print(colSums(selected_features))
	num_features <- as.numeric(colSums(selected_features))
	# print("num_features:")
	# print(num_features)
	# print("length(num_features):")
	# print(length(num_features))
	problem_lambda_indices <- 1:ncol(selected_features)
	problem_lambda_indices <- problem_lambda_indices[num_features==0]
	# print("problem_lambda_indices:")
	# print(problem_lambda_indices)
	problem_lambdas <- unique(selected_features_lambdas[problem_lambda_indices])
	# print("problem_lambdas:")
	# print(problem_lambdas)
	# Identify all sets of feature sets corresponding to these lambdas and
	# remove them
	# print("problem lambdas:")
	# print(problem_lambdas)
	if(length(problem_lambdas) == 1){
		# Indices of sets to remove
		sets_to_keep <- selected_features_lambdas!=as.numeric(problem_lambdas)
		# print("sets_to_keep:")
		# print((1:ncol(selected_features))[sets_to_keep])
		# print("sets to remove:")
		# print((1:ncol(selected_features))[!sets_to_keep])
		# print("sets to remove:")
		# print(sets_to_remove)
		selected_features <- selected_features[, sets_to_keep]
		selected_features_lambdas <- selected_features_lambdas[sets_to_keep]
		num_features <- as.numeric(colSums(selected_features))
	} else if(length(problem_lambdas) > 1){
		sets_to_remove <- selected_features_lambdas %in% problem_lambdas
		sets_to_keep <- !sets_to_remove
		# for(i in 1:length(problem_lambdas)){
		# 	sets_to_remove <- c(sets_to_remove, (0:B)*nlambda +
		# 		problem_lambdas[i])
		# }
		# sets_to_remove <- sets_to_remove[sets_to_remove <= B*nlambda]
		# sets_to_remove <- sets_to_remove[sets_to_remove > 0]
		# print("sets to remove:")
		# print(sets_to_remove)
		selected_features <- selected_features[, sets_to_keep]
		selected_features_lambdas <- selected_features_lambdas[sets_to_keep]
		num_features <- as.numeric(colSums(selected_features))
		# if(!is.matrix(selected_features) || !is.numeric(selected_features)){
	 #    	print("Error2: selected_features is not a matrix.")
	 #    }
	 #    if(ncol(selected_features) != length(selected_features_lambdas)){
	 #    	print("ERROR: ncol(selected_features) != length(selected_features_lambdas)")
	 #    }
	}
	# print("as.numeric(colSums(selected_features)) (should be fixed now):")
	# print(as.numeric(colSums(selected_features)))
	# print("length(num_features):")
	# # print(length(num_features))
	# print("ncol(selected_features):")
	# print(ncol(selected_features))
	# print("length(selected_features_lambdas):")
	# print(length(selected_features_lambdas))
	return(list(selected_features, selected_features_lambdas, num_features))
}

calculateCanonicalCorrelations <- function(X.dat, selected_features,
	selected_features_k, digest_table, num_feat_sets, num_lambdas){
	require(ccaPP)
	# For each set of selected features, calculate all of the the maximum
	# canonical correlations between the selected feature set and the
	# remaining selected feature sets. Take the median of these canonical
	# correlations as our metric for closeness to the center.
	# Inputs:
	# 
	# X.dat: n x p design matrix
	#
	# selected_features: a p x num_feat_sets matrix of selected features.
	# Each column contains a set of selected features, with a 1 in row i of
	# column j if feature i was selected in iteration j, and a 0 otherwise.
	#
	# selected_features_k: a p x num_feats[k] matrix containing only the
	# columns of selected_features corresponding to the current lambda.
	#
	# digest_table: hashmap. Each time a canonical 
	# correlation between two feature sets is calculated, an entry will
	# be added to this hashmap. The key of this entry will be the
	# digest number for the pair of feature sets, and the value will
	# be the canonical correlation between those two feature sets.
	#
	# num_feat_sets: the total number of feature sets under consideration (B for
	# each lambda that did not have an empty feature set on any resample.)
	#
	# num_lambdas:  Number of lambdas under consideration (number of lambda
	# values fitted in lasso paths)
	# 
	# 
	# Outputs:
	#
	# cancors_k: a num_feat_sets by num_feat_sets_k matrix; each of the
	# num_feat_sets_k columns will contain the num_feat_sets maximum canonical 
	# correlations the ith feature space (for this lambda) has with each of the
	# num_feat_sets feature spaces.
	#
	# digest_table: the digest_table with added entries for any new canonical
	# correlations calculated

	# Number of feature sets under consideration (corresponding to the current
	# lambdas)
	num_feat_sets_k <- ncol(selected_features_k)
	if(is.null(num_feat_sets_k)){
		num_feat_sets_k <- 1
	}

	# if(is.vector(selected_features_k) && !is.list(selected_features_k)){
	# 	num_feat_sets_k <- 1
	# }
	# else{
	# 	num_feat_sets_k <- ncol(selected_features_k)
	# 	}
	# if(!is.numeric(num_feat_sets_k)){
	# 	print("ERROR: !is.numeric(num_feat_sets_k")
	# 	print(num_feat_sets_k)
	# 	print(selected_features_k)
	# }
	# if(!is.numeric(num_feat_sets)){
	# 	print("ERROR: !is.numeric(num_feat_sets")
	# 	print(num_feat_sets)
	# }

	cancors_k <- matrix(0, num_feat_sets, num_feat_sets_k)

	# cat("starting the double for loop...", fill=TRUE)
	for(i in 1:num_feat_sets_k){
		# print(paste("i = ", i, " out of ", B))
		for(j in 1:num_feat_sets){
			if(num_feat_sets_k > 1){
				feat_set_1 <- selected_features_k[, i]
			} else{
				feat_set_1 <- selected_features_k
			}
			
			feat_set_2 <- selected_features[, j]
			# Get digest_key for these two feature sets
			digest_key_1 <- digest(list(feat_set_1, feat_set_2))
			# Check if the canonical correlation between these
			# two selected feature sets has already been calculated;
			# if so, return it
			# if (digest_table$has_key(digest_key_1)){
			# 	cancors_k[j, i] <- digest_table$find(digest_key_1)
			# }
			value <- digest_table$find(digest_key_1)
			if (!is.na(value)){
				cancors_k[j, i] <- value
			}
			# else if (digest_key_2 %in% digest_table$Digests){
			# 	cancors_k[j, i] <-
			# 		digest_table$CCs[digest_table$Digests==digest_key_2]
			# } 
			else{
				# If the canonical correlation between these two
				# feature sets has not already been calculated, calculate
				# it now and add it to the digest_table
		 		cancor <- ccaPP::maxCorGrid(X.dat[, feat_set_1],
		 			X.dat[, feat_set_2], method="pearson")$cor
		 		cancors_k[j, i] <- cancor
		 		if(is.na(cancor)){
		 			print("Error: cancor is NA. wtf")
		 			print("feat_set_1:")
		 			# print(feat_set_1)
		 			print((1:nrow(selected_features))[feat_set_1])
		 			print("feat_set_2:")
		 			# print(feat_set_2)
		 			print((1:nrow(selected_features))[feat_set_2])
		 		}
		 		# digest_table[nrow(digest_table) + 1, ] <- list(
		 		# 	digest_key_1, cancors_k[j, i])
		 		digest_key_2 <- digest(list(feat_set_2, feat_set_1))
		 		digest_table$insert(c(digest_key_1, digest_key_2),
		 			rep(cancor, 2))
		 		# print(paste(i, ", ", j))
			}
		}
	}
	return(list(cancors_k, digest_table))
}

getSelectedFeaturePathCanCor <- function(X.dat, selected_features,
	selected_features_lambdas, B){
	# Inputs: 
	#
	# X.dat: n x p design matrix
	#
	# selected_features: a p x num_feat_sets matrix of selected features.
	# Each column contains a set of selected features, with a 1 in row i of
	# column j if feature i was selected in iteration j, and a 0 otherwise.
	#
	# selected_features_lambdas: an ncol(selected_features)-vector that
	# keeps track of the lambda corresponding to each selected feature set in
	# selected_features. The ith entry is the lambda corresponding to the
	# selected feature set in the ith column of selected_features.
	#
	# B: number of resamples used. Should also be the number of selected
	# feature sets at each lambda.
	#
	# Output: 
	# 
	# selected: a list of selected features of length num_lambdas.
	# Each selected feature set will be a vector.
	#
	# digest_table: a hashmap containing digests of the pairs of feature sets
	# as keys and maximum canonical correlations between the two feature sets
	# as values.

	selected <- list()
	# Digest table: hashmap. Each time a canonical 
	# correlation between two feature sets is calculated, an entry will
	# be added to this hashmap. The key of this entry will be the
	# digest number for the pair of feature sets, and the value will
	# be the canonical correlation between those two feature sets.
	digest_table <- hashmap(character(1), numeric(1))
	# Call the total number of feature sets under consideration
	# num_feat_sets
	# num_feat_sets: the total number of feature sets under consideration (B for
	# each lambda that did not have an empty feature set on any resample.)
	num_feat_sets <- ncol(selected_features)
	# Unique lambdas under consideration
	lambdas <- sort(unique(selected_features_lambdas))
	# num_lambdas:  Number of lambdas under consideration (number of lambda
	# values fitted in lasso paths)
	num_lambdas <- length(lambdas)
	# colnames(digest_table) <- c("Digests", "CCs")
	# For each lambda, do the following:
	for(k in 1:num_lambdas){
		# print(paste(k, " out of ", num_lambdas))
		# print(Sys.time() - t1)

		# Create a new matrix that only includes the columns of 
		# selected_features corresponding to this lambda
		sets_to_consider_indices <- selected_features_lambdas==lambdas[k]
		if(length(sets_to_consider_indices) == 0){
			print("ERROR: length(sets_to_consider_indices) == 0")
		}
		selected_features_k <- selected_features[, sets_to_consider_indices]
		
		# For each set of selected features, calculate all of the the maximum
		# canonical correlations between the selected feature set and the
		# remaining selected feature sets. 
			
		cancors_k_results <- calculateCanonicalCorrelations(X.dat=X.dat,
			selected_features=selected_features,
			selected_features_k=selected_features_k, digest_table=digest_table,
			num_feat_sets=num_feat_sets, num_lambdas=num_lambdas)
		cancors_k <- cancors_k_results[[1]]
		# print("cancors_k:")
		# print(cancors_k)
		digest_table <- cancors_k_results[[2]]
		rm(cancors_k_results)
		# Now cancors_k is a num_feat_sets by B matrix; each of the B columns
		# contains the num_feat_sets maximum canonical correlations the ith
		# feature space (for this lambda) has with each of the num_feat_sets
		# feature spaces.

		# Take the median of these canonical
		# correlations as our metric for closeness to the center.
		# median_cancors_k is a B-dimensional vector that contains the median
		# maximum canonical correlation each of the B feature spaces under 
		# consideration has with the sum_i nlambda_i feature spaces.
		median_cancors_k <- apply(cancors_k, 2, median)
		print("median_cancors_k:")
		print(median_cancors_k)
		# The feature set with the largest median canonical correlation
		# is the one we want (the best one for this lambda).
		stability_selected_indices_k <- median_cancors_k==max(median_cancors_k)
		# If there is a tie between highest median canonical correlation,
		# arbitratily choose the first feature space.
		stability_selected_index_k <-
			min(which(stability_selected_indices_k))
		# print(stability_selected_index_k)
		if(!is.null(ncol(selected_features_k))){
			stability_selected_k  <- as.numeric(which(selected_features_k[,
				stability_selected_index_k] !=0))
		} else{
			stability_selected_k  <- as.numeric(which(selected_features_k[
				stability_selected_index_k] !=0))
		}

		
		# Add these selections to the list of selections
		selected[[k]] <- stability_selected_k
	}
	return(list(selected, digest_table))
}


lambda.max.cancor.ss <- new_method(name = "lambda.max.cancor.ss",
	label = "L_M_CC_SS",
	method = function(model, draw, B, nlambda) {
		t1 <- Sys.time()
		n <- nrow(model$x)
		p <- ncol(model$x)
		# Fit B resamples of lasso on lambda paths of size nlambda
		bs_results <- bootstrapSamplingIndicatorsQMax(X.dat=model$x,
			y=draw, B=B, nlambda=nlambda)
		selected_features <- bs_results[[1]]
		n_feature_sets <- bs_results[[2]]
		selected_features_lambdas <- bs_results[[3]]
		rm(bs_results)
		# Now n_feature_sets contains the number of selected feature sets for
		# each of the B lasso fits, and selected_features contains
		# B*sum(n_feature_sets) columns of selected features. 
		# Remove any lambdas that contained at least one column of all zeros.
		zerolambdas_results <- removeZeroLambdas(selected_features=
			selected_features,
			selected_features_lambdas=selected_features_lambdas, 
			nlambda=nlambda, B=B)
		selected_features <- zerolambdas_results[[1]]
		# print("selected_features")
		# print(selected_features)
		# print(dim(selected_features))
		# print(p)
		# print(B*nlambda)
		selected_features_lambdas <- zerolambdas_results[[2]]
		# print("selected_features_lambdas:")
		# print(selected_features_lambdas)
		num_features <- zerolambdas_results[[3]]
		rm(zerolambdas_results)
		# num_features <- as.numeric(colSums(selected_features))
		
		# We will return a list of selected features, one at each lambda.
		selected_results <- getSelectedFeaturePathCanCor(X.dat=model$x,
			selected_features=selected_features,
			selected_features_lambdas=selected_features_lambdas, B=B)
		selected <- selected_results[[1]]
		# print("selected:")
		# print(selected)
		# Now selected is a list of selected features of length num_lambdas.
		digest_table <- selected_results[[2]]
		rm(selected_results)
		# Lastly, to make size eval work, create a matrix of betas which are
		# just 1 if the coefficient is nonzero and zero otherwise.
		beta <- matrix(0, p, length(selected))
		for(i in 1:length(selected)){
			beta[selected[[i]], i] <- 1
		}
		# print("Complete")
		print("Finished a sample")
		print(Sys.time() - t1)
		# print("Final length of digest_table (divided by 2):")
		# print((digest_table$size()-1)/2)
		# print("Selected features:")
		# print(selected)
		return(list(selected=selected, beta=beta))
	},
	settings=list(B=20, nlambda=25)
)
