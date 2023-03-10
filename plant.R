# setwd("/Users/gregfaletto/Documents/GitHub/css-sims")

library(cssr)
library(hdf5r)
library(dplyr)
library(glmnet)
library(Metrics)
library(ggplot2)
library(parallel)
library(doParallel)
library(ggcorrplot)

library(gridExtra)
library(cowplot)
library(simulator)

cl <- parallel::makeForkCluster(nnodes=detectCores())
doParallel::registerDoParallel(cl)

# https://dept.stat.lsa.umich.edu/~jerrick/courses/stat701/notes/parallel.html
# https://privefl.github.io/blog/a-guide-to-parallelism-in-r/
# http://www.john-ros.com/Rcourse/parallel.html

sim_dir <- getwd()

# load data?
load_data <- FALSE

# Run new study, or load study that has been previously run?
run_new_study <- FALSE

# Training set proportion

selec_prop <- 0.4
train_prop <- 0.4

# Correlation cutoff for clustering
cor_cutoff <- 0.5

# Number of draws to take
n_draws <- 1000
# n_draws <- 200

# Number of SNPs to use in data set
n_snps <- 1000

# Largest model size to use in evaluations
p_max <- 60

# Largest model size to use in plots
p_max_plots <- 60

stopifnot(p_max_plots <= p_max)

# Coarseness for stability plots (how many model sizes to include in one point?)
# coarseness <- round(p_max/20)
coarseness <- 5

# Minimum number of observations to take an average for MSE or calculate
# stability metric
MIN_COUNT <- round(n_draws*.01*coarseness)

# Verbose printing in loops?
verbose <- FALSE

# Directory where simFunctions.R is stored
wd <- getwd()
dir_funcs <- paste(sim_dir, "Helper Functions", sep="/")
# dir_resp <- "/Users/gregfaletto/Google Drive/Data Science/LaTeX/Paper skeleton 2020/Real data examples/GWAS/AraGWAS/database/12"
# dir_hdf5 <- "/Users/gregfaletto/Google Drive/Data Science/Python/gwas"
dir_resp <- "/Users/gregfaletto/My Drive/Data Science/LaTeX/Paper skeleton 2020/Real data examples/GWAS/AraGWAS/database/12"
dir_hdf5 <- "/Users/gregfaletto/My Drive/Data Science/Python/gwas"

# SNPs file name
snps_file <- "snps50000.csv"

############# Load Functions #################

setwd(dir_funcs)
source(file="simFunctions.R")
source("toy_ex_slide_funcs.R")
source("model_functions.R")
source("method_functions.R")
source("eval_functions.R")
setwd(wd)

if(load_data){

	t_data <- Sys.time()

	# Load SNP data
	setwd(dir_hdf5)

	df <- H5File$new("4.hdf5", mode="r")
	df$ls(recursive=TRUE)

	accessions <- as.integer(df[["accessions"]][])
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

	stopifnot(ncol(snps) <= length(positions))
	colnames(snps) <- paste("x", positions[1:ncol(snps)], sep="")

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

	# • position was not in genetic map. **** doesn't apply

	print("Done! Number of SNPs:")

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

	stopifnot(nrow(pheno) == length(unique(pheno$accession_id)))

	data_plant <- as.data.frame(snps)

	data_plant$accession_id <- as.integer(rownames(data_plant))
	data_plant <- left_join(data_plant, pheno, by="accession_id")
	data_plant <- data_plant[, colnames(data_plant) != "accession_id"]

	if(any(is.na(data_plant))){
		stop("any(is.na(data_plant))")
	}

	rownames(snps) <- paste("accession", rownames(snps), sep="")

	print("Total time to load and process data:")
	print(Sys.time() - t_data)
}

setwd(wd)

t_max <- Sys.time()

if(run_new_study){

	t1 <- Sys.time()

	response_name <- "log_FT10"

	if(!(response_name %in% colnames(data_plant))){
		stopifnot("FT10" %in% colnames(data_plant))
		# Change response to log
		data_plant$FT10 <- log(data_plant$FT10)
		colnames(data_plant)[colnames(data_plant) == "FT10"] <- response_name
	}

	stopifnot(response_name %in% colnames(data_plant))
	
	# Selection, training and test sets

	n <- nrow(data_plant)
	n_selec <- round(selec_prop*n)
	n_train <- round(train_prop*n)
	n_test <- n - n_train - n_selec

	set.seed(341)

	plant_sim <- new_simulation("plant_sim", "GWAS Simulation Study") |>
        generate_model(make_model=plant_model, n_selec=n_selec, n_train=n_train,
        	n_test=n_test, cor_cutoff=cor_cutoff, response_name=response_name,
        	max_model_size=p_max) |> simulate_from_model(nsim = n_draws) |>
        run_method(c(SS_SS_cssr_plant # Stability selection (as proposed by Shah and
            # Samworth 2012)
            , SS_SS_cssr_elnet_plant # Stability selection with elastic net as
            # base procedure
            , SS_CSS_sparse_cssr_plant # Sparse cluster stability selection
            , SS_CSS_weighted_cssr_plant # Weighted averaged cluster stability
            # selection
            , SS_CSS_avg_cssr_plant # Simple averaged cluster stability
            # selection
            , clusRepLasso_cssr_plant # Cluster representative lasso
            , protolasso_cssr_plant # Protolasso
            , lasso_random_plant # Lasso
            , elastic_net_plant # Elastic Net
        ))

    save_simulation(plant_sim)

    plant_sim <- plant_sim |> evaluate(list(cssr_mse_plant))

    save_simulation(plant_sim)

    print("Total time for simulation:")
    print(Sys.time() - t1)
    print("Time per simulation:")
    print((Sys.time() - t1)/n_draws)

} else{
	# Load needed files
	plant_sim <- load_simulation("plant_sim")
}

parallel::stopCluster(cl)

### Generate figures

results <- genPlotDfPlant(plant_sim, coarseness=coarseness)

results_df <- results$results_df

n_methods <- results$n_methods

e_df <- results$eval_df

# Table of means of MSEs

sum_res_avg <- df_plant_app(e_df, max_model_size=p_max,
    methods=c("SS_CSS_avg_cssr_plant", "SS_CSS_weighted_cssr_plant",
        "SS_SS_cssr_plant", "clusRepLasso_cssr_plant", "protolasso_cssr_plant",
        "lasso_random_plant"), coarseness=coarseness)

print("")
print("Means of MSES:")
print("")

stargazer(sum_res_avg, summary=FALSE)

# Table of means of stability

stab_df <- plantStabSummary(results_df, max_model_size=p_max,
	num_sims=n_draws, methods=nameMap(c("SS_CSS_avg_cssr_plant",
		"SS_CSS_weighted_cssr_plant",
        "SS_SS_cssr_plant", "clusRepLasso_cssr_plant", "protolasso_cssr_plant",
        "lasso_random_plant")), coarseness=coarseness)

print("")
print("Means of NSB Stability:")
print("")

stargazer(stab_df, summary=FALSE)

### Figure 4 (previously Figure 5) 
# Omit uncompetitive methods for visual clarity

fig_4_left <- createLossesPlot3(results_df[!(results_df$Method %in%
    nameMap(c("SS_CSS_sparse_cssr_plant", "elastic_net_plant",
    	"SS_SS_cssr_elnet_plant"))), ], n_methods - 3, plot_errors=FALSE,
	max_model_size=p_max, log_mse=TRUE, break_by=2*coarseness)

fig_4_mid <- createNSBStabPlot2(results_df[!(results_df$Method %in%
    nameMap(c("SS_CSS_sparse_cssr_plant", "elastic_net_plant",
    	"SS_SS_cssr_elnet_plant"))), ], plot_errors=FALSE, break_by=2*coarseness)

fig_4_right <- createStabMSEPlot2(results_df[!(results_df$Method %in%
    nameMap(c("SS_CSS_sparse_cssr_plant", "elastic_net_plant",
    	"SS_SS_cssr_elnet_plant"))), ], n_methods - 3, plot_errors=FALSE,
	log_mse=TRUE)

# 2. Save the legend
#+++++++++++++++++++++++
legend <- get_legend(fig_4_left + theme(legend.direction="horizontal"))

# 3. Remove the legend from the box plot
#+++++++++++++++++++++++
fig_4_left <- fig_4_left + theme(legend.position="none")

fig_4_mid <- fig_4_mid + theme(legend.position="none")

fig_4_right <- fig_4_right + theme(legend.position="none")

# 4. Arrange ggplot2 graphs with a specific width

fig_4 <- grid.arrange(fig_4_left, fig_4_mid, fig_4_right, legend, ncol=3,
    nrow = 2, layout_matrix = rbind(c(1, 2, 3), c(4, 4, 4)),
    widths = c(1.8, 1.8, 1.8), heights = c(2.5, 0.2))

fig_4 <- cowplot::ggdraw(fig_4) +
    theme(plot.background = element_rect(fill="white", color = NA))

print(fig_4)

saveFigure2(subdir="figures", plot=fig_4, size="large",
	filename="fig_real_data.pdf")

### Versions of Figure 4 plots with all methods (for supplement)

fig_4_supp_left <- createLossesPlot3(results_df, n_methods,
	plot_errors=FALSE, max_model_size=p_max, log_mse=TRUE, break_by=coarseness)

saveFigure2(subdir="figures", plot=fig_4_supp_left, size="xmlarge",
    filename="real_data_mse_supp.pdf")

fig_4_supp_mid <- createNSBStabPlot2(results_df, plot_errors=FALSE,
	break_by=coarseness)

saveFigure2(subdir="figures", plot=fig_4_supp_mid, size="xmlarge",
    filename="real_data_stab_supp.pdf")

fig_4_supp_right <- createStabMSEPlot2(results_df, n_methods, plot_errors=FALSE,
	log_mse=TRUE)

saveFigure2(subdir="figures", plot=fig_4_supp_right, size="xmlarge",
    filename="real_data_mse_stab_supp.pdf")

print("Done with simulations! Total time:")
print(Sys.time() - t_max)
print("Time per simulation:")
print((Sys.time() - t_max)/n_draws)








