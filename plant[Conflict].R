# setwd("/Users/gregfaletto/Google Drive/Data Science/LaTeX/Paper skeleton 2020/Real data examples/GWAS/AraGWAS/GENOTYPES")

# dev.off()
rm(list=ls())

library(hdf5r)
library(dplyr)
library(glmnet)
library(Metrics)
library(ggplot2)

# https://stat.ethz.ch/pipermail/r-help/2016-December/443661.html
Sys.setenv("MC_CORES"=8L)

library(doMC)
library(parallel)
library(doParallel)
# library(rhdf5)
# library(rhdf5filters)

# dev.off()
rm(list=ls())

registerDoMC(cores = detectCores())
registerDoParallel()
options(mc.cores = parallel::detectCores()-1)

# Training set proportion

selec_prop <- 0.4
train_prop <- 0.4
if(selec_prop + train_prop >= 1){
	stop("selec_prop + train_prop >= 1")
}
# n_train <- 200

# Number of draws to take
n_draws <- 10

# Largest model size to use in evaluations
p_max <- 100

# Verbose prrinting in loops?
verbose <- FALSE

# Directory where simFunctions.R is stored
dir_funcs <- "/Users/gregfaletto/Google Drive/Data Science/LaTeX/Paper skeleton 2020/Real data examples"

dir_funcs2 <- "/Users/gregfaletto/Google Drive/Data Science/LaTeX/Generalized Stability Selection Presentation"


# Directory where stabsel_copy.R is storeed
dir_stabsel <- "/Users/gregfaletto/Google Drive/Data Science/R/USC/Stability Selection/toy_example/New method"

# Methods
methods <- c("lasso"
	, "ss_ss" # Shah and Samwoth stability selection
	, "ss_gss" # Generalized stability selection (Shah and Samworth version)
	, "ss_gss_avg" # Generalized stability selection (Shah and Samworth version)
	# with weighted averaging
	, "ss_gss_avg_unwt" # Generalized stability selection (Shah and Samworth version)
	# with unweighted (simple) averaging
	)

n_methods <- length(methods)

dir_main <- getwd()

dir_pheno <- "/Users/gregfaletto/Google Drive/Data Science/LaTeX/Paper skeleton 2020/Real data examples/GWAS/AraGWAS/database/12"

dir_hdf5 <- "/Users/gregfaletto/Google Drive/Data Science/Python/gwas"


# SNPs file name

snps_file <- "snps50000.csv"







############# Load Functions #################

setwd(dir_stabsel)
source(file="stabsel_copy.R")

setwd(dir_funcs)
source(file="simFunctions.R")

setwd(dir_funcs2)
source("toy_ex_slide_funcs.R")

setwd(dir_main)

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

setwd(dir_main)

# Load responses (phenotypes)

setwd(dir_pheno)

pheno <- read.csv("study_12_values.csv", header=TRUE)

setwd(dir_main)

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
setwd(dir_main)
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

colnames(snps) <- positions[1:ncol(snps)]

# SNPs satisfying any of the following conditions were removed:
# • minor allele frequency (MAF) < 1%,
freqs <- colSums(snps)/nrow(snps)
print(paste("Removing", sum(freqs < 0.01 | freqs > 0.99),
	"SNPs with a minor allele frequency (MAF) < 1%."))
snps <- snps[, freqs >= 0.01 & freqs <= 0.99]
# • Hardy–Weinberg equilibrium test p-value < 0.01%,
# https://www.rdocumentation.org/packages/genetics/versions/1.3.8.1.2/topics/HWE.test
# https://cran.r-project.org/web/packages/HardyWeinberg/

# • missing > 5% of values (values were considered missing if the BRLMM score 
# was > 0.5), ***already taken care of--imputed data

# • position was listed as same as another SNP ***** doesn't apply

# • position was not in genetic map. **** look into this? probably fine




print("Done!")

print("Number of SNPs:")

p_snps <- ncol(snps)

print(p_snps)

if(p_snps > 1000){
	print("Reducing to 1,000 SNPs...")
	snps <- snps[, 1:1000]
	p_snps <- 1000
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

# List of lists of Matrices to store selected sets in, as binary vectors (in
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

	# Random draw
	selec_indices <- sample(1:n, n_selec)
	train_indices <- sample(setdiff(1:n, selec_indices), n_train)
	test_indices <- setdiff(1:n, c(selec_indices, train_indices))

	data_selec <- data[selec_indices, ]
	data_train <- data[train_indices, ]
	data_test <- data[test_indices, ]

	# # Create matrices

	# data_mat_train <- as.matrix(data[train_indices, colnames(data) != "FT10"])
	# data_mat_test <- as.matrix(data[test_indices, colnames(data) != "FT10"])

	# y_train <- data[train_indices, colnames(data) == "FT10"]
	# y_test <- data[test_indices, colnames(data) == "FT10"]

	# Get clusters and put in R matrix

	R <- getR(snps_mat=snps[c(selec_indices, train_indices), ], verbose=verbose)

	# Get selected snps from various methods using training set

	for(k in 1:n_methods){
		# Get ranked list of features to select
		ranked_results_k <- getRankedFeatures(data_selec, methods[k], j_choices,
			R)
		selections_k <- ranked_results_k[[1]]
		j_choices <- ranked_results_k[[2]]
		
		if(length(j_choices) == 0){
            print(paste("length(j_choices) == 0 on iteration", i))
            next
        }

        print(paste("Finished getting selected sets for method", methods[k]))
        print(Sys.time() - t0)

		for(j in j_choices){
			# print("selections_k[[j]]:")
			# print(selections_k[[j]])
			if(length(selections_k[[j]]) == 0){
				stop("length(selections_k[[j]]) == 0")
			}
			if(!all(selections_k[[j]] %in% colnames(snps))){
				stop("!all(selections_k[[j]] %in% colnames(snps))")
			}
			if(sum(colnames(snps) %in% selections_k[[j]]) == 0) {
				stop("sum(colnames(snps) %in% selections_k[[j]]]) == 0")
			}
			sel_mats[[k]][[j]][i, colnames(snps) %in% selections_k[[j]]] <- 1
		}

		# Keep track of these j_choices
        j_choices_mat[j_choices + (k-1)*p_max, i] <- TRUE
        
        # Calculate MSEs

        for(j in j_choices){
        	if(sum(sel_mats[[k]][[j]][i, ] == 1) == 0){
        		stop("sum(sel_mats[[k]][[j]][i, ] == 1) == 0)")
        	}
        	selected_kj <- colnames(snps)[sel_mats[[k]][[j]][i, ] == 1]
        	if(methods[k] %in% c("ss_gss_avg", "ss_gss_avg_unwt")){
                losses_mat[(k - 1)*p_max + j, i] <- getLossPlant(x_train,
                    y_train, x_test, y_test, selected_kj,
                    average=TRUE, to_avg=to_avg_list_list[[k]][[j]],
                    avg_feats=avg_feats_list_list[[k]][[j]],
                    weights=weights_list_list[[k]][[j]])
            } else{
            	losses_mat[(k - 1)*p_max + j, i] <- getLossPlant(data_train,
                	data_test, selected_kj)
            }
        }
	}

	print(paste("Time for iteration", i, ":"))
	print(Sys.time() - t0)
}

# MSE Plot

mse_plot <- createLossesPlot3(losses_mat, methods, j_choices_mat, legend=TRUE)

print(mse_plot)


# Finally, calculate stability metrics and create plot of those.
stab_mets_Nogueira <- matrix(as.numeric(NA), p_max, n_methods)
stab_mets_Nogueira_conf_lower <- matrix(as.numeric(NA), p_max,
    n_methods)
stab_mets_Nogueira_conf_upper <- matrix(as.numeric(NA), p_max,
    n_methods)
for(j in 1:p_max){
    for(k in 1:n_methods){
        stat_results <- calcNogueiraStab(sel_mats[[k]][[j]], n_draws,
            calc_errors=TRUE)
        stab_mets_Nogueira[j, k] <- stat_results[1]
        stab_mets_Nogueira_conf_lower[j, k] <- stat_results[2]
        stab_mets_Nogueira_conf_upper[j, k] <- stat_results[3]
    }
}

colnames(stab_mets_Nogueira) <- methods

# print("stab_mets_Nogueira:")
# print(stab_mets_Nogueira)

plot_stability_Nogueira_ints <- createNogueiraStabPlot3(stab_mets_Nogueira,
    plot_errors = TRUE, lowers=stab_mets_Nogueira_conf_lower,
    uppers=stab_mets_Nogueira_conf_upper)

print(plot_stability_Nogueira_ints)

plot_stability_Nogueira <- createNogueiraStabPlot3(stab_mets_Nogueira,
    plot_errors = FALSE)

print(plot_stability_Nogueira)




