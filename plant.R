# setwd("/Users/gregfaletto/Documents/GitHub/css-sims")
# setwd("/Users/gregfaletto/Google Drive/Data Science/LaTeX/Paper skeleton 2020/Real data examples/GWAS/AraGWAS/GENOTYPES")

library(hdf5r)
library(dplyr)
library(glmnet)
library(Metrics)
library(ggplot2)
# Sys.setenv("MC_CORES"=8L)
# https://dept.stat.lsa.umich.edu/~jerrick/courses/stat701/notes/parallel.html
# library(doMC)
library(parallel)
library(doParallel)
library(ggcorrplot)

library(gridExtra)
library(cowplot)
# library(rhdf5)
# library(rhdf5filters)

# dev.off()
# rm(list=ls())

# registerDoMC(cores = detectCores())
# registerDoParallel()
# options(mc.cores = parallel::detectCores()-1)

cl <- parallel::makeForkCluster(nnodes=detectCores())
doParallel::registerDoParallel(cl)

# https://dept.stat.lsa.umich.edu/~jerrick/courses/stat701/notes/parallel.html
# https://privefl.github.io/blog/a-guide-to-parallelism-in-r/
# http://www.john-ros.com/Rcourse/parallel.html

sim_dir <- getwd()

# load data?
load_data <- FALSE

# Conduct exploratory analysis?
exp_analysis <- FALSE

# Run new study, or load study that has been previously run?
run_new_study <- FALSE

# Generate plots?
gen_plots <- TRUE

# Training set proportion

selec_prop <- 0.4
train_prop <- 0.4
if(selec_prop + train_prop >= 1){
	stop("selec_prop + train_prop >= 1")
}
# n_train <- 200

# Correlation cutoff for clustering (for now, pick either 0.1 or 0.5)
# cor_cutoff <- 0.1
cor_cutoff <- 0.5

folder <- paste("210523/Real Data Example (cor_cutoff = ", cor_cutoff ,
	", small sample size)", 
	# ")",
	sep="")


folder_main <- "/Users/gregfaletto/Dropbox/Subspace Stability Selection/sims_6"
folder_dir <- file.path(folder_main, folder)
dir.create(folder_dir, showWarnings = FALSE, recursive = TRUE)

# Number of draws to take
n_draws <- 100
# n_draws <- 75
# n_draws <- 2

# Number of SNPs to use in data set
n_snps <- 1000

# Largest model size to use in evaluations
p_max <- round(n_snps/10)

# Largest model size to use in plots
p_max_plots <- round(n_snps/10)
# p_max_plots <- 60

stopifnot(p_max_plots <= p_max)

# Coarseness for stability plots (how many model sizes to include in one point?)
coarseness <- round(p_max/20)

coarseness_plots <- round(p_max_plots/20)

# Verbose printing in loops?
verbose <- FALSE


# Directory where simFunctions.R is stored
wd <- getwd()
dir_funcs <- paste(sim_dir, "Helper Functions", sep="/")
dir_resp <- "/Users/gregfaletto/Google Drive/Data Science/LaTeX/Paper skeleton 2020/Real data examples/GWAS/AraGWAS/database/12"
dir_hdf5 <- "/Users/gregfaletto/Google Drive/Data Science/Python/gwas"
# dir_funcs <- "/Users/gregfaletto/Google Drive/Data Science/LaTeX/Paper skeleton 2020/Real data examples"

# dir_funcs2 <- "/Users/gregfaletto/Google Drive/Data Science/LaTeX/Generalized Stability Selection Presentation"


# # Directory where stabsel_copy.R is storeed
# dir_stabsel <- "/Users/gregfaletto/Google Drive/Data Science/R/USC/Stability Selection/toy_example/New method"

# Methods
methods <- c("lasso"
	, "lasso_proto" # Protolasso
	, "SS_SS" # Shah and Samwoth stability selection
	, "SS_GSS" # Sparse generalized stability selection (Shah and Samworth version)
	, "SS_GSS_avg" # Generalized stability selection (Shah and Samworth version)
	# with weighted averaging
	, "SS_GSS_avg_unwt" # Generalized stability selection (Shah and Samworth version)
	# with unweighted (simple) averaging
	, "BRVZ_avg_unwt" # Cluster representative lasso
	)

n_methods <- length(methods)

# Averaging methods
averaging_methods <- c("SS_GSS_avg", "SS_GSS_avg_unwt", "BRVZ_avg_unwt")



# dir_pheno <- "/Users/gregfaletto/Google Drive/Data Science/LaTeX/Paper skeleton 2020/Real data examples/GWAS/AraGWAS/database/12"




# SNPs file name

snps_file <- "snps50000.csv"







############# Load Functions #################

# setwd(dir_stabsel)
# source(file="stabsel_copy.R")

setwd(dir_funcs)
source(file="simFunctions.R")

# setwd(dir_funcs2)
source("toy_ex_slide_funcs.R")

setwd(wd)


if(load_data){

	# Load SNP data

	# data <- hdf5r::readDataSet("4.hdf5")


	# h5ls(file="4.hdf5")
	# data <- rhdf5::h5read(file="4.hdf5", name="snps")

	# lzf filter issue: https://support.bioconductor.org/p/130731/

	setwd(dir_hdf5)

	df <- H5File$new("4.hdf5", mode="r")

	df$ls(recursive=TRUE)

	# str(df[["accessions"]][])

	accessions <- as.integer(df[["accessions"]][])

	# str(df[["snps"]][1:2, ])

	# # workaround I need: https://github.com/satijalab/seurat/issues/2977#issuecomment-629351772

	# str(df[["positions"]][])

	positions <- df[["positions"]][]

	df$close_all()

	setwd(wd)

	# Load responses (phenotypes)

	setwd(dir_resp)

	pheno <- read.csv("study_12_values.csv", header=TRUE)

	setwd(wd)

	print("Number of unique accessions in response data:")

	print(length(unique(pheno$accession_id)))

	print("Number of these that are in SNP data:")

	print(sum(unique(pheno$accession_id) %in% accessions))

	# Acceessions to get

	acc_inds <- which(accessions %in% pheno$accession_id)


	# Read snps file
	print("Loading SNPs data...")
	t0 <- Sys.time()
	setwd(dir_hdf5)
	snps <- t(read.csv(snps_file, header=FALSE))
	setwd(wd)
	print("Done! Total time to load:")
	print(Sys.time() - t0)

	if(nrow(snps) != length(accessions)){
		stop("nrow(snps) != length(accessions)")
	}

	rownames(snps) <- accessions

	# Eliminate accesssions not in SNP data
	print("Keeping only accessions that have responses (phenotypes) available...")
	snps <- snps[as.integer(rownames(snps)) %in% pheno$accession_id, ]
	print("Done!")

	if(ncol(snps) > length(positions)){
		stop("ncol(snps) > length(positions)")
	}

	colnames(snps) <- paste("x", positions[1:ncol(snps)], sep="")
	# print("colnames(snps)[1:10]:")
	# print(colnames(snps)[1:10])

	# print(unique(apply(snps, 2, function(x) length(unique(x)))))

	# unique(apply(snps, 2, function(x) length(unique(x))))

	# SNPs satisfying any of the following conditions were removed:
	# • minor allele frequency (MAF) < 1%,
	freqs <- colSums(snps)/nrow(snps)
	print(paste("Removing", sum(freqs < 0.01 | freqs > 0.99),
		"SNPs with a minor allele frequency (MAF) < 1%."))
	snps <- snps[, freqs >= 0.01 & freqs <= 0.99]
	# Double-check
	sds <- apply(snps, 2, sd)
	snps <- snps[, sds > 0]
	# • Hardy–Weinberg equilibrium test p-value < 0.01%,
	# https://www.rdocumentation.org/packages/genetics/versions/1.3.8.1.2/topics/HWE.test
	# https://cran.r-project.org/web/packages/HardyWeinberg/
	# https://en.wikipedia.org/wiki/Hardy%E2%80%93Weinberg_principle#Significance_tests_for_deviation

	# • missing > 5% of values (values were considered missing if the BRLMM score 
	# was > 0.5), ***already taken care of--imputed data

	# • position was listed as same as another SNP ***** doesn't apply

	# • position was not in genetic map. **** look into this? probably fine




	print("Done!")

	print("Number of SNPs:")

	p_snps <- ncol(snps)

	print(p_snps)

	if(p_snps > n_snps){
		print(paste("Reducing to", n_snps, "SNPs..."))
		snps <- snps[, 1:n_snps]
		p_snps <- n_snps
		print("Done!")
	}



	# Get response (use FT10--no NAs, and also duplicate entry for accession 8424
	# has same value in both replicates)

	pheno <- pheno[pheno$accession_id %in% as.integer(rownames(snps)),
		c("accession_id", "FT10")]

	pheno <- pheno[!duplicated(pheno), ]

	if(nrow(pheno) != length(unique(pheno$accession_id))){
		stop("nrow(pheno) != length(unique(pheno$accession_id))")
	}

	data <- as.data.frame(snps)

	data$accession_id <- as.integer(rownames(data))

	data <- left_join(data, pheno, by="accession_id")

	data <- data[, colnames(data) != "accession_id"]

	if(any(is.na(data))){
		stop("any(is.na(data))")
	}

	rownames(snps) <- paste("accession", rownames(snps), sep="")
}

if(exp_analysis){
	
	# 	Look at correlation matrix—say first 500 features
	cor_mat <- cor(snps)

    print("Done! Time to calculate:")
    print(Sys.time() - t0)

    # Plot correlations

    no_names_cor_mat <- cor_mat

    colnames(no_names_cor_mat) <- character()
    rownames(no_names_cor_mat) <- character()

    print(ggcorrplot(no_names_cor_mat[1:50, 1:50]))

    cors_df <- as.data.frame(no_names_cor_mat[lower.tri(no_names_cor_mat,
    	diag = FALSE)])

    colnames(cors_df) <- "Correlations"

    cor_plot <- ggplot(cors_df, aes(x=Correlations)) + geom_histogram()

    print(cor_plot)

    # rm(cor_plot)

    print(summary(cors_df))

    print("Proportion of absolute correlations greater than 0.5:")
    print(sum(abs(cors_df$Correlations) > 0.5)/length(cors_df$Correlations))

    print("Proportion of absolute correlations greater than 0.9:")
    print(sum(abs(cors_df$Correlations) > 0.9)/length(cors_df$Correlations))

    print("str(snps):")
    print(str(snps))
   
    


    # print(ggcorrplot(cor_mat))

    # print(sort(abs(cor_mat[, colnames(cor_mat)=="gdp_growth"]),
    #   decreasing=TRUE))

    dist <- as.dist(1 - abs(cor_mat))
    h <- hclust(dist)
    # rm(cor_mat)


    # Look at dendrogram
    plot(h)
    # View(data)

    # we clustered the SNPs using the estimated correlations as a similarity
    # measure with a single-linkage cutoff of 0.5, and settle for discovering 
    # important SNP clusters.
    ct <-  cutree(h, h=cor_cutoff)
    # ct
    # table(ct)

    print("Number of clusters of each size:")

    print(paste("Number of clusters estimated in full data set (cutoff =",
    	cor_cutoff, "):"))

    print(length(unique(ct)))

    print("Number of SNPs:")

    print(ncol(snps))

    print("Average number of SNPs per cluster:")

    print(ncol(snps)/length(unique(ct)))

    print("Number of clusters of each size:")

    print(table(table(ct)))


	# 	Form hierarchical clustering based on Pearson correlation like Candes (but more stringent clustering)



	# 	Histogram of response—should we log transform? Summary statistics
	response_plot <- ggplot(data, aes(x=FT10)) + geom_histogram()
	print(response_plot)

	log_response_plot <- ggplot(data, aes(x=log(FT10))) + geom_histogram()
	print(log_response_plot)
	# 	Summary statistics of X matrix—histogram of frequencies at which different SNPs are present
	# What proportion of samples have a given mutation
	freqs <- colMeans(snps)
	freqs <- data.frame(freqs)
	colnames(freqs) <- "Frequencies"
	freq_hist <- ggplot(freqs, aes(x=Frequencies)) + geom_histogram()
	print(freq_hist)
	# •	
	# 	Histogram of number of mutations in a given plant
}

if(run_new_study){

	response_name <- "log_FT10"

	if(!(response_name %in% colnames(data))){
		if(!("FT10" %in% colnames(data))){
		stop("!(FT10 %in% colnames(data))")
	}
		# Change response to log
		data$FT10 <- log(data$FT10)
		colnames(data)[colnames(data) == "FT10"] <- response_name
	}

	if(!(response_name %in% colnames(data))){
		stop("!(response_name %in% colnames(data))")
	}
	
	# Selection, training and test sets

	n <- nrow(data)
	n_selec <- round(selec_prop*n)
	n_train <- round(train_prop*n)
	n_test <- n - n_train - n_selec

	# # Matrix of MSE values
	# mse_mat <- matrix(NA, n_draws, n_methods)
	# colnames(mse_mat) <- methods

	# Matrix to store MSE results in

	losses_mat <- matrix(0, nrow=n_methods*p_max, ncol=n_draws)
	losses_mat_rows <- character()
	for(i in 1:n_methods){
	    loss_name_i <- paste(methods[i], "loss")
	    losses_mat_rows <- c(losses_mat_rows, paste(loss_name_i, 1:p_max))
	}
	if(length(losses_mat_rows) != nrow(losses_mat)){
	    stop("length(losses_mat_rows) != nrow(losses_mat)")
	}
	rownames(losses_mat) <- losses_mat_rows

	# List of lists of matrices to store selected sets in, as binary vectors (in
	# order to calculate stability metrics)

	sel_mats <- list()
	for(k in 1:n_methods){
	    list_k <- list()
	    for(j in 1:p_max){
	    	mat_kj <- matrix(0, n_draws, p_snps)
	    	colnames(mat_kj) <- colnames(snps)
	        list_k[[j]] <- mat_kj
	    }
	    sel_mats[[k]] <- list_k
	    rm(list_k)
	    rm(mat_kj)
	}

	# List of lists of matrices to store selected sets in, as binary vectors (in
	# order to calculate stability metrics)--for purpose of stability metric
	# (i.e., for averaging methods, includes all selected features in selected
	# cluster)

	sel_mats_stab <- list()
	for(k in 1:n_methods){
	    list_k <- list()
	    for(j in 1:p_max){
	    	mat_kj <- matrix(0, n_draws, p_snps)
	    	colnames(mat_kj) <- colnames(snps)
	        list_k[[j]] <- mat_kj
	    }
	    sel_mats_stab[[k]] <- list_k
	    rm(list_k)
	    rm(mat_kj)
	}



	# method_indices <- integer(n_methods)

	# for(k in 1:n_methods){
	#     method_indices[k] <- identifyMethodIndex(methods[k],
	#         output_p_prime)
	# }

	# Matrix to store which model sizes existed for each sim
	j_choices_mat <- matrix(FALSE, nrow=n_methods*p_max, ncol=n_draws)

	# Set seed for reproducibility
	set.seed(341)

	for(i in 1:n_draws){

		t0 <- Sys.time()

		# Possible model sizes
		j_choices <- 1:p_max

		# If using averaging method, need lists for features to average,
	    # clusters to average, and weights
	    if(any(averaging_methods %in% methods)){
	        to_avg_list_list <- list()
	        avg_feats_list_list <- list()
	        weights_list_list <- list()
	    }

		# Random draw
		selec_indices <- sample(1:n, n_selec)
		train_indices <- sample(setdiff(1:n, selec_indices), n_train)
		test_indices <- setdiff(1:n, c(selec_indices, train_indices))

		data_selec <- data[selec_indices, ]
		data_train <- data[train_indices, ]
		data_test <- data[test_indices, ]

		snps_i <- as.matrix(data_selec[, colnames(data_selec) != response_name])

		# SNPs satisfying any of the following conditions were removed:
		# • minor allele frequency (MAF) < 1%,
		freqs_i <- colSums(snps_i)/nrow(snps_i)
		inds_to_remove_i <- which(freqs_i < 0.01 | freqs_i > 0.99)
		# if(length(inds_to_remove_i) > 0){
		# 	print(paste("Removing", length(inds_to_remove_i),
		# 		"SNPs with a minor allele frequency (MAF) < 1%."))
		# 	snp_inds <- setdiff(1:ncol(snps), inds_to_remove_i)
		# 	var_names <- colnames(snps)[snp_inds]
		# 	inds_to_keep_i <- setdiff(1:ncol(data_selec), inds_to_remove_i)
		# 	data_selec <- data_selec[, inds_to_keep_i]
		# 	data_train <- data_train[, inds_to_keep_i]
		# 	data_test <- data_test[, inds_to_keep_i]
		# } else{
			var_names <- colnames(snps)
			snp_inds <- 1:ncol(snps)
		# }
		
		# # Double-check
		# sds <- apply(snps, 2, sd)
		# snps <- snps[, sds > 0]

		# # Create matrices

		# data_mat_train <- as.matrix(data[train_indices, colnames(data) != "FT10"])
		# data_mat_test <- as.matrix(data[test_indices, colnames(data) != "FT10"])

		# y_train <- data[train_indices, colnames(data) == "FT10"]
		# y_test <- data[test_indices, colnames(data) == "FT10"]

		# Get clusters and put in R matrix

		R <- getR(snps_mat=snps[c(selec_indices, train_indices), snp_inds],
			cor_cutoff=cor_cutoff, verbose=verbose)

		if(ncol(R) != length(var_names)){
			stop("ncol(R) != length(var_names)")
		}

		# Get selected snps from various methods using training set

		for(k in 1:n_methods){
			t1 <- Sys.time()
			# Get ranked list of features to select
			ranked_results_k <- getRankedFeatures(data_selec, methods[k],
				j_choices, var_names=var_names, response_name=response_name,
				R=R)
			if(all(is.null(ranked_results_k[[1]]))){
	            stop("all(is.null(ranked_results_k[[1]]))")
	        }
			selections_k <- ranked_results_k$selected_sets
			j_choices_k <- ranked_results_k$j_choices
			
			if(length(j_choices_k) == 0){
	            print(paste("length(j_choices_k) == 0 on iteration", i))
	            next
	        }

	        if(methods[k] %in% averaging_methods){
	            to_avg_list_list[[k]] <- ranked_results_k$to_avg_list
	            avg_feats_list_list[[k]] <- ranked_results_k$avg_feats_list
	            weights_list_list[[k]] <- ranked_results_k$weights_list

	            # Outputs:
	            #
	            # to_avg_list_list: A list of length n_meths_to_eval. Each
	            # element is a list to_avg_list whose elements are logical
	            # vectors of length j. For method k, if feature l in 1:j is in a 
	            # cluster, the lth entry of the jth vector of the kth
	            # to_avg_list in to_avg_list_list will be TRUE. 
	            #
	            # avg_feats_list_list: A list of length n_meths_to_eval. Each
	            # element is a list avg_feats_list whose elements are lists of
	            # length j. For method k, if feature l in 
	            # 1:j is in a cluster, the lth entry of the jth vector of the kth
	            # avg_feats_list in avg_feats_list_list
	            # will contain a character vector containing the names of the 
	            # features from the cluster of which feature l is a member.
	            #
	            # weights_list_list: A list of length n_meths_to_eval. Each
	            # element is a list weights_list whose elements are lists of
	            # length j. For method k, if feature l in 
	            # 1:j is in a cluster, the lth entry of the jth vector of the
	            # kth weights_list in weights_list_list will
	            # be the weights to use (in the same order as the jth entry of
	            # the kth avg_feats_list in avg_feats_list_list).

	        }

	        # print(paste("Finished getting selected sets for method",
	        # 	nameMap(methods[k])))
	        # print("Total time for feature selection and extraction:")
	        # print(Sys.time() - t1)

			for(j in j_choices_k){
				# print("selections_k[[j]]:")
				# print(selections_k[[j]])
				if(length(selections_k[[j]]) == 0){
					stop("length(selections_k[[j]]) == 0")
				}
				if(!all(selections_k[[j]] %in% colnames(snps))){
					inds_bad <- !(selections_k[[j]] %in% colnames(snps))
					print("bad selections:")
					print(selections_k[[j]][inds_bad])
					print("good selections:")
					print(selections_k[[j]][!inds_bad])
					stop("!all(selections_k[[j]] %in% colnames(snps))")
				}
				if(sum(colnames(snps) %in% selections_k[[j]]) == 0) {
					stop("sum(colnames(snps) %in% selections_k[[j]]]) == 0")
				}
				sel_mats[[k]][[j]][i, colnames(snps) %in% selections_k[[j]]] <- 1
				sel_mats_stab[[k]][[j]][i, colnames(snps) %in% selections_k[[j]]] <- 1
				if(methods[k] %in% averaging_methods){
					# Are any of the selected features cluster members?

                    # If none of the selected features are cluster members,
                    # the cluster representation is the same as the non-
                    # cluster representation

                    if(any(to_avg_list_list[[k]][[j]])){
                        # Identify the cluster members
                        clust_feats <- which(to_avg_list_list[[k]][[j]])
                        for(l in clust_feats){
                            selected_feat_inds_k_j_l <- avg_feats_list_list[[k]][[j]][[l]]
                            sel_mats_stab[[k]][[j]][i, colnames(snps) %in%
                            	selected_feat_inds_k_j_l] <- 1
                        }
                    } 
				}
			}

			# Keep track of these j_choices
	        j_choices_mat[j_choices_k + (k-1)*p_max, i] <- TRUE
	        
	        # Calculate MSEs

	        for(j in j_choices_k){
	        	if(sum(sel_mats[[k]][[j]][i, ] == 1) == 0){
	        		stop("sum(sel_mats[[k]][[j]][i, ] == 1) == 0)")
	        	}
	        	selected_kj <- colnames(snps)[sel_mats[[k]][[j]][i, ] == 1]
	        	if(methods[k] %in% averaging_methods){
	                losses_mat[(k - 1)*p_max + j, i] <- getLossPlant(data_train,
	                	data_test, selected_kj, var_names=var_names,
	                	response_name=response_name,
	                    average=TRUE, to_avg=to_avg_list_list[[k]][[j]],
	                    avg_feats=avg_feats_list_list[[k]][[j]],
	                    weights=weights_list_list[[k]][[j]])
	            } else{
	            	losses_mat[(k - 1)*p_max + j, i] <- getLossPlant(data_train,
	                	data_test, selected_kj, var_names=var_names,
	                	response_name=response_name)
	            }
	        }
		}

		print(paste("Time for iteration", i, ":"))
		print(Sys.time() - t0)
	}

	# Save data

	save(losses_mat, file=paste("losses_mat_", cor_cutoff, ".Rdata", sep=""))

	save(methods, file=paste("methods_", cor_cutoff, ".Rdata", sep=""))

	save(j_choices_mat, file=paste("j_choices_mat_", cor_cutoff, ".Rdata",
		sep=""))

	save(sel_mats, file=paste("sel_mats_", cor_cutoff, ".Rdata", sep=""))

	save(sel_mats_stab, file=paste("sel_mats_stab_", cor_cutoff, ".Rdata", sep=""))

	save(n_draws, file=paste("n_draws_", cor_cutoff, ".Rdata", sep=""))

	save(p_max, file=paste("p_max_", cor_cutoff, ".Rdata", sep=""))
} else{
	# Load needed files

	load(paste("losses_mat_", cor_cutoff, ".Rdata", sep=""))

	load(paste("methods_", cor_cutoff, ".Rdata", sep=""))

	load(paste("j_choices_mat_", cor_cutoff, ".Rdata",
		sep=""))

	load(paste("sel_mats_", cor_cutoff, ".Rdata", sep=""))

	load(paste("sel_mats_stab_", cor_cutoff, ".Rdata", sep=""))

	load(paste("n_draws_", cor_cutoff, ".Rdata", sep=""))

	load(paste("p_max_", cor_cutoff, ".Rdata", sep=""))
}

if(gen_plots){

	# MSE Plot

	mse_plot <- createLossesPlot3(losses_mat, methods, j_choices_mat,
		legend=TRUE, coarseness=coarseness_plots, p_max_plots=p_max_plots)

	print(mse_plot)

	saveFigure("MSE vs Num Fitted Coefs", mse_plot, size="small",
		filename="MSE vs Num Fitted Coefs.pdf")

	if(p_max/coarseness != round(p_max/coarseness)){
		stop("p_max/coarseness != round(p_max/coarseness)")
	}


	# Finally, calculate stability metrics and create plot of those.
	stab_mets_NSB <- matrix(as.numeric(NA), p_max/coarseness, n_methods)
	# stab_mets_NSB_conf_lower <- matrix(as.numeric(NA), p_max/coarseness,
	# 	n_methods)
	# stab_mets_NSB_conf_upper <- matrix(as.numeric(NA), p_max/coarseness,
	# 	n_methods)
	for(j in 1:(p_max/coarseness)){
	    for(k in 1:n_methods){
	    	# Assemble sel_mats_stab I want
	    	if(coarseness > 1){
				j_indices <- ((j-1)*coarseness + 1):(j*coarseness)
				sel_mat_kj <- sel_mats_stab[[k]][[j_indices[1]]]
				for(j_prime in 2:length(j_indices)){
					sel_mat_kj <- rbind(sel_mat_kj,
						sel_mats_stab[[k]][[j_indices[j_prime]]])
				}
				# print("dim(sel_mat_kj):")
				# print(dim(sel_mat_kj))
				# print("sel_mat_kj:")
				# print(sel_mat_kj)
			} else{
				sel_mat_kj <- sel_mats_stab[[k]][[j]]
			}
	        stab_mets_NSB[j, k] <- calcNSBStab(sel_mat_kj, n_draws,
	        	calc_errors=FALSE, coarseness=coarseness)
	        # stab_mets_NSB[j, k] <- stat_results[1]
	        # stab_mets_NSB_conf_lower[j, k] <- stat_results[2]
	        # stab_mets_NSB_conf_upper[j, k] <- stat_results[3]
	    }
	}

	colnames(stab_mets_NSB) <- methods

	# print("stab_mets_NSB:")
	# print(stab_mets_NSB)

	# plot_stability_NSB_ints <- createNSBStabPlot3(stab_mets_NSB,
	#     coarseness=coarseness, plot_errors = TRUE,
	#     lowers=stab_mets_NSB_conf_lower, uppers=stab_mets_NSB_conf_upper)

	# print(plot_stability_NSB_ints)

	plot_stability_NSB <- createNSBStabPlot3(stab_mets_NSB,
	    coarseness=coarseness, plot_errors = FALSE)

	print(plot_stability_NSB)

	saveFigure("Stability vs Num Fitted Coefficients/Original Stabilty Metric",
        plot_stability_NSB, size="small",
        filename="Stability vs Num Fitted Coefficients.pdf")

	stab_mse_plot <- createStabMSEPlot3(losses=losses_mat,
		stab_mets=stab_mets_NSB, names=methods, j_choices_mat=j_choices_mat,
		p_max_plots=p_max_plots, legend=TRUE, coarseness=coarseness,
		plot_errors=FALSE
		# , names_to_omit=c("lasso_proto", "BRVZ_avg_unwt")
		)

	print(stab_mse_plot)

	saveFigure("Stability vs MSE/Original Stabilty Metric", stab_mse_plot,
		size="medium", filename="Stability vs MSE.pdf")

	saveFigure(subdir="Slides", plot=stab_mse_plot, size="slide",
        filename="real_data_stab_vs_mse.pdf")











	# Side-by-side plot

    # 2. Save the legend
    #+++++++++++++++++++++++
    legend <- get_legend(mse_plot + theme(legend.direction="horizontal"))
    # 3. Remove the legend from the box plot
    #+++++++++++++++++++++++
    mse_plot_no_legend <- mse_plot + theme(legend.position="none")
    plot_stability_NSB_no_legend <- plot_stability_NSB +
        theme(legend.position="none")
    blankPlot <- ggplot() + geom_blank(aes(1,1)) + cowplot::theme_nothing()
    # 4. Arrange ggplot2 graphs with a specific width

    # side_plot <- grid.arrange(mse_plot_no_legend,
    #     plot_stability_NSB_no_legend, legend
    #     , ncol=3
    #     , widths=c(0.4, 0.4, 0.2)
    #     )

    # side_plot <- cowplot::ggdraw(side_plot) +
    #     theme(plot.background = element_rect(fill="white", color = NA))

    # print(side_plot)

    # saveFigure("MSE and Stability vs Num Fitted Coefficients/Original Stabilty Metric/Confidence Intervals",
    #     side_plot, size="large")



	# 3 figure Side-by-side plot

    # 3. Remove the legend from the box plot
    #+++++++++++++++++++++++
    stab_mse_plot_no_legend <- stab_mse_plot + theme(legend.position="none")

    # 4. Arrange ggplot2 graphs with a specific width

    side_plot <- grid.arrange(mse_plot_no_legend, plot_stability_NSB_no_legend,
    	stab_mse_plot_no_legend, legend, ncol=3, nrow = 2,
    	layout_matrix = rbind(c(1, 2, 3), c(4, 4, 4)),
    	widths = c(1.8, 1.8, 1.8), heights = c(2.5, 0.2))

    side_plot <- cowplot::ggdraw(side_plot) +
        theme(plot.background = element_rect(fill="white", color = NA))

    print(side_plot)

    saveFigure("MSE and Stability vs Num Fitted Coefficients/Original Stabilty Metric",
        side_plot, size="xlarge", filename="3_side_plot.pdf")

}

parallel::stopCluster(cl)


