## @knitr models

# make_my_model <- function(n, prob) {
#   new_model(name = "contaminated-normal",
#             label = sprintf("Contaminated normal (n = %s, prob = %s)", n, prob),
#             params = list(n = n, mu = 2, prob = prob),
#             simulate = function(n, mu, prob, nsim) {
#               # this function must return a list of length nsim
#               contam <- runif(n * nsim) < prob
#               x <- matrix(rep(NA, n * nsim), n, nsim)
#               x[contam] <- rexp(sum(contam))
#               x[!contam] <- rnorm(sum(!contam))
#               x <- mu + x # true mean is mu
#               return(split(x, col(x))) # make each col its own list element
#             })
# }







make_coefficients <- function(p, k, beta_min, beta_max, nblocks, sig_blocks,
    block_size){
    # true_coefs <- runif(k, min=beta_min, max=beta_max)
    beta <- numeric(p)
    # identify indices of first coefficient in each significant block (these
    # betas will be nonzero)
    blocked_dgp_vars <- ((0:(sig_blocks-1))*block_size+1)
    # make these coefficients in each block nonzero
    beta[blocked_dgp_vars] <- runif(sig_blocks, min=beta_min, max=beta_max)
    # identify remaining coefficients in blocks (which ought to be set to 0)
    insig_blocked_vars <- setdiff(1:(block_size*nblocks), blocked_dgp_vars)
    # find significant unblocked variables (if applicable) and fill in 
    # coefficients
    if(k > sig_blocks){
        sig_unblocked_vars <- (nblocks*block_size + 1):(nblocks*block_size + k -
            sig_blocks)
        beta[sig_unblocked_vars] <- runif(k - sig_blocks, min=beta_min,
            max=beta_max)
    } else {
        sig_unblocked_vars <- NA
    }
    return(list(beta, blocked_dgp_vars, sig_unblocked_vars, insig_blocked_vars))
}


make_sparse_blocked_linear_model <- function(n, p, k, beta_min,
    beta_max, nblocks, sig_blocks, block_size, rho, var, snr) {
    # x <- matrix(rnorm(n * p), n, p)
    Sigma <- make_covariance_matrix(p, nblocks, block_size, rho, var)
    x <- mvrnorm(n=n, mu=rep(0, p), Sigma=Sigma)
    coefs <- make_coefficients(p, k, beta_min, beta_max, nblocks, sig_blocks,
        block_size)
    beta <- coefs[[1]]
    # blocked_dgp_vars <- coefs[2]
    # sig_unblocked_vars <- coefs[3]
    mu <- as.numeric(x %*% beta)
    sd <- sqrt(sum(mu^2) / (n * snr)) # taking snr = ||mu||^2 / (n * sigma^2)
    my_model <- new_model(  name = "sblm", 
                label = sprintf("Linear model with correlated blocks (n = %s,
                    p = %s, k = %s, nblocks = %s, sig_blocks = %s, rho = %s)",
                    n, p, k, nblocks, sig_blocks, rho),
                params = list(x = x, beta = beta, mu = mu, sd = sd, n = n,
                    p = p, k = k, beta_min = beta_min, beta_max = beta_max,
                    nblocks = nblocks, sig_blocks = sig_blocks,
                    block_size = block_size, rho = rho, var = var,
                    blocked_dgp_vars = coefs[[2]],
                    sig_unblocked_vars = coefs[[3]],
                    insig_blocked_vars = coefs[[4]]),
                simulate = function(mu, sd, nsim) {
                    y <- mu + sd * matrix(rnorm(nsim * n), n, nsim)
                    return(split(y, col(y))) #make each col its own list element
                }
    )
    # print("finished generating model")
    return(my_model)
}

make_coefficients2 <- function(p, k, beta_low, beta_high, nblocks, sig_blocks,
    block_size){
    # Same as make_coefficients, but sets blocked coefficients to beta_low and
    # unblocked coefficients to beta_high
    beta <- numeric(p)
    # identify indices of first coefficient in each significant block (these
    # betas will be nonzero)
    blocked_dgp_vars <- ((0:(sig_blocks-1))*block_size+1)
    # make these coefficients in each block nonzero
    beta[blocked_dgp_vars] <- beta_low
    # identify remaining coefficients in blocks (which ought to be set to 0)
    insig_blocked_vars <- setdiff(1:(block_size*nblocks), blocked_dgp_vars)
    # find significant unblocked variables (if applicable) and fill in 
    # coefficients
    if(k > sig_blocks){
        sig_unblocked_vars <- (nblocks*block_size + 1):(nblocks*block_size + k -
            sig_blocks)
        beta[sig_unblocked_vars] <- rep(beta_high, k - sig_blocks)
    } else {
        sig_unblocked_vars <- NA
    }
    return(list(beta, blocked_dgp_vars, sig_unblocked_vars, insig_blocked_vars))
}

make_coefficients4 <- function(p, k_unblocked=0, beta_low, beta_high, nblocks,
    sig_blocks, block_size){
    # sets blocked coefficients to beta_high and
    # unblocked coefficients to beta_low
    beta <- numeric(p)
    # identify indices of first coefficient in each significant block (these
    # betas will be nonzero)
    blocked_dgp_vars <- ((0:(sig_blocks-1))*block_size+1)
    # make these coefficients in each block nonzero
    beta[blocked_dgp_vars] <- beta_high
    # identify remaining coefficients in blocks (which ought to be set to 0)
    insig_blocked_vars <- setdiff(1:(block_size*nblocks-1), blocked_dgp_vars)
    # find significant unblocked variables (if applicable) and fill in 
    # coefficients
    if(k_unblocked > 0){
        sig_unblocked_vars <- (nblocks*block_size + 1):
        (nblocks*block_size + k_unblocked)
        beta[sig_unblocked_vars] <- rep(beta_low, k_unblocked)
    } else {
        sig_unblocked_vars <- NA
    }
    # print("sig_unblocked_vars:")
    # print(sig_unblocked_vars)
    return(list(beta, blocked_dgp_vars, sig_unblocked_vars, insig_blocked_vars))
}

make_coefficients4_ranking <- function(p, k_unblocked=0, beta_low_min,
    beta_low_max, beta_high, nblocks, sig_blocks, block_size){
    # sets blocked coefficients to beta_high and
    # unblocked coefficients to range from beta_low_max to beta_low_min
    beta <- numeric(p)
    # identify indices of first coefficient in each significant block (these
    # betas will be nonzero)
    blocked_dgp_vars <- ((0:(sig_blocks-1))*block_size+1)
    # make these coefficients in each block nonzero
    beta[blocked_dgp_vars] <- beta_high
    # identify remaining coefficients in blocks (which ought to be set to 0)
    insig_blocked_vars <- setdiff(1:(block_size*nblocks-1), blocked_dgp_vars)
    # Range of weak signal coefficients
    beta_lows <- seq(beta_low_max, beta_low_min, length.out=k_unblocked)
    # find significant unblocked variables (if applicable) and fill in 
    # coefficients
    if(k_unblocked > 0){
        sig_unblocked_vars <- (nblocks*block_size + 1):
        (nblocks*block_size + k_unblocked)
        beta[sig_unblocked_vars] <- beta_lows
    } else {
        sig_unblocked_vars <- NA
    }
    # print("sig_unblocked_vars:")
    # print(sig_unblocked_vars)
    return(list(beta, blocked_dgp_vars, sig_unblocked_vars, insig_blocked_vars))
}



make_sparse_blocked_linear_model2 <- function(n, p, k, beta_low,
    beta_high, nblocks, sig_blocks, block_size, rho, var, snr) {
    # Same as make_sparse_blocked_linear_model, but sets blocked coefficients to
    # beta_low and unblocked coefficients to beta_high, and also removes the
    # true signal feature from x after generating the response y.
    Sigma <- make_covariance_matrix(p, nblocks, block_size, rho, var)
    x <- mvrnorm(n=n, mu=rep(0, p), Sigma=Sigma)
    coefs <- make_coefficients2(p, k, beta_low, beta_high, nblocks, sig_blocks,
        block_size)
    beta <- coefs[[1]]
    # blocked_dgp_vars <- coefs[2]
    # sig_unblocked_vars <- coefs[3]
    mu <- as.numeric(x %*% beta)
    # Remove true blocked signal feature from x now that I've generated mu
    x <- x[, 2:ncol(x)]
    sd <- sqrt(sum(mu^2) / (n * snr)) # taking snr = ||mu||^2 / (n * sigma^2)
    my_model <- new_model(  name = "sblm2", 
                label = sprintf("Linear model2 with correlated blocks (n = %s,
                    p = %s, k = %s, nblocks = %s, sig_blocks = %s, rho = %s)",
                    n, p, k, nblocks, sig_blocks, rho),
                params = list(x = x, beta = beta, mu = mu, sd = sd, n = n,
                    p = p, k = k, beta_low = beta_low, beta_high = beta_high,
                    nblocks = nblocks, sig_blocks = sig_blocks,
                    block_size = block_size, rho = rho, var = var,
                    blocked_dgp_vars = coefs[[2]],
                    sig_unblocked_vars = coefs[[3]],
                    insig_blocked_vars = coefs[[4]]),
                simulate = function(mu, sd, nsim) {
                    y <- mu + sd * matrix(rnorm(nsim * n), n, nsim)
                    return(split(y, col(y))) #make each col its own list element
                }
    )
    # print("finished generating model")
    return(my_model)
}



make_coefficients3 <- function(p, k, beta_low, beta_high, nblocks, sig_blocks,
    block_size){
    # Same as make_coefficients, but sets blocked coefficients to beta_low and
    # unblocked coefficients to beta_high
    beta <- numeric(p)
    # identify indices of first coefficient in each significant block (these
    # betas will be nonzero)
    if(sig_blocks > 0){
        blocked_dgp_vars <- ((0:(sig_blocks-1))*block_size+1)
        # make these coefficients in each block nonzero
        beta[blocked_dgp_vars] <- beta_high
    } else{
        blocked_dgp_vars <- numeric()
    }
    # identify remaining coefficients in blocks (which ought to be set to 0)
    if(sig_blocks > 0){
        insig_blocked_vars <- setdiff(1:(block_size*nblocks), blocked_dgp_vars)
    } else{
        insig_blocked_vars <- 1:(block_size*nblocks)
    }
    # find significant unblocked variables (if applicable) and fill in 
    # coefficients
    if(k > sig_blocks){
        sig_unblocked_vars <- (nblocks*block_size + 1):(nblocks*block_size + k -
            sig_blocks)
        beta[sig_unblocked_vars] <- rep(beta_low, k - sig_blocks)
    } else {
        sig_unblocked_vars <- NA
    }
    return(list(beta, blocked_dgp_vars, sig_unblocked_vars, insig_blocked_vars))
}

make_sparse_blocked_linear_model3 <- function(n, p, k, beta_low,
    beta_high, nblocks, sig_blocks, block_size, rho, var, snr) {
    # Same as make_sparse_blocked_linear_model, but sets blocked coefficients to
    # beta_high and unblocked coefficients to beta_low. Does not remove the
    # true signal feature from x after generating the response y.

    # x <- matrix(rnorm(n * p), n, p)
    Sigma <- make_covariance_matrix(p, nblocks, block_size, rho, var)
    x <- mvrnorm(n=n, mu=rep(0, p), Sigma=Sigma)
    coefs <- make_coefficients3(p, k, beta_low, beta_high, nblocks, sig_blocks,
        block_size)
    beta <- coefs[[1]]
    # blocked_dgp_vars <- coefs[2]
    # sig_unblocked_vars <- coefs[3]
    mu <- as.numeric(x %*% beta)
    sd <- sqrt(sum(mu^2) / (n * snr)) # taking snr = ||mu||^2 / (n * sigma^2)
    my_model <- new_model(  name = "sblm", 
                label = sprintf("Linear model with correlated blocks (n = %s,
                    p = %s, k = %s, nblocks = %s, sig_blocks = %s, rho = %s)",
                    n, p, k, nblocks, sig_blocks, rho),
                params = list(x = x, beta = beta, mu = mu, sd = sd, n = n,
                    p = p, k = k, beta_low = beta_low, beta_high = beta_high,
                    nblocks = nblocks, sig_blocks = sig_blocks,
                    block_size = block_size, rho = rho, var = var,
                    Sigma = Sigma,
                    blocked_dgp_vars = coefs[[2]],
                    sig_unblocked_vars = coefs[[3]],
                    insig_blocked_vars = coefs[[4]]),
                simulate = function(mu, sd, nsim) {
                    y <- mu + sd * matrix(rnorm(nsim * n), n, nsim)
                    return(split(y, col(y))) #make each col its own list element
                }
    )
    # print("finished generating model")
    return(my_model)
}

gen_mu_x_sd4 <- function(n, p, k_unblocked, beta_low, beta_high, nblocks=1,
    sig_blocks=1, block_size, rho, var, snr=NA, sigma_eps_sq=NA){
    if(is.na(snr) & is.na(sigma_eps_sq)){
        stop("Must specify one of snr or sigma_eps_sq")
    }
    Sigma <- make_covariance_matrix(p + sig_blocks, nblocks, block_size +
        sig_blocks, rho, var)
    # print(Sigma)
    x <- mvrnorm(n=n, mu=rep(0, p + sig_blocks), Sigma=Sigma)
    coefs <- make_coefficients4(p + sig_blocks, k_unblocked, beta_low,
        beta_high, nblocks, sig_blocks, block_size + 1)
    beta <- coefs[[1]]
    # blocked_dgp_vars <- coefs[2]
    # sig_unblocked_vars <- coefs[3]
    mu <- as.numeric(x %*% beta)
    # print(beta)
    # Remove true blocked signal feature from x now that I've generated mu
    # (for now only have code written for one block)
    if(nblocks==1 & sig_blocks ==1){
        z <- x[, 1]
        x <- x[, 2:ncol(x)]
    } else{
        stop("!(nblocks==1 & sig_blocks ==1)")
    }
    # If SNR is null, use sigma_eps_sq
    if(!is.na(sigma_eps_sq)){
        sd <- sqrt(sigma_eps_sq)
    }else{
        sd <- sqrt(sum(mu^2) / (n * snr)) # taking snr = ||mu||^2 /(n * sigma^2)
    }
    # Take out first row of Sigma that no longer needs to be there (again, 
    # assuming one block)
    if(nblocks==1 & sig_blocks ==1){
        Sigma <- Sigma[2:nrow(Sigma), 2:ncol(Sigma)]
    }
    return(list(mu=mu, x=x, sd=sd, Sigma=Sigma, beta=beta, z=z, coefs=coefs))
}


make_sparse_blocked_linear_model4 <- function(n, p, k_unblocked, beta_low,
    beta_high, nblocks=1, sig_blocks=1, block_size, rho, var, snr=NA,
    sigma_eps_sq=NA) {
    # Same as make_sparse_blocked_linear_model, but sets blocked coefficients to
    # beta_high and unblocked coefficients to beta_low, and also removes the
    # true signal feature from x after generating the response y.
    if(is.na(snr) & is.na(sigma_eps_sq)){
        stop("Must specify one of snr or sigma_eps_sq")
    }
    mu_x_sd <- gen_mu_x_sd4(n, p, k_unblocked, beta_low, beta_high, nblocks,
        sig_blocks, block_size, rho, var, snr, sigma_eps_sq)
    mu <- mu_x_sd[[1]]
    x <- mu_x_sd[[2]]
    sd <- mu_x_sd[[3]]
    Sigma <- mu_x_sd[[4]]
    beta <- mu_x_sd[[5]]
    coefs <- mu_x_sd[[7]]
    rm(mu_x_sd)
    
    my_model <- new_model(  name = "sblm2", 
                label = sprintf("Linear model2 with correlated blocks (n = %s,
                    p = %s, k_unblocked = %s, nblocks = %s, sig_blocks = %s, rho = %s)",
                    n, p, k_unblocked, nblocks, sig_blocks, rho),
                params = list(x = x, beta = beta, mu = mu, sd = sd, n = n,
                    p = p, k_unblocked = k_unblocked, beta_low = beta_low,
                    beta_high = beta_high,
                    nblocks = nblocks, sig_blocks = sig_blocks,
                    block_size = block_size, rho = rho, var = var,
                    Sigma = Sigma,
                    blocked_dgp_vars = coefs[[2]],
                    sig_unblocked_vars = coefs[[3]],
                    insig_blocked_vars = coefs[[4]]),
                simulate = function(mu, sd, nsim) {
                    y <- mu + sd * matrix(rnorm(nsim * n), n, nsim)
                    return(split(y, col(y))) #make each col its own list element
                }
    )
    # print("finished generating model")
    return(my_model)
}


random_simulate_func <- function(n, p, k_unblocked, beta_low, beta_high,
    nblocks=1, sig_blocks=1, block_size, rho, var, snr=NA, sigma_eps_sq=NA,
    Sigma, beta, nsim){
    if(is.na(snr) & is.na(sigma_eps_sq)){
        stop("Must specify one of snr or sigma_eps_sq")
    }

    # List we'll return: will be length nsim, and every element will be a
    # named list with elements X and y.
    ret_list <- list()
    stopifnot(nrow(Sigma) == p + sig_blocks & ncol(Sigma) == p + sig_blocks)
    if(!all(dim(Sigma) == c(p + sig_blocks, p + sig_blocks))){
        print("dim(Sigma):")
        print(dim(Sigma))
        print("p:")
        print(p)
        print("sig_blocks:")
        print(sig_blocks)
        stop("nrow(Sigma) != p + sig_blocks | ncol(Sigma) != p + sig_blocks")
    }
    for(i in 1:nsim){
        # print(Sigma)
        mu_set <- rep(0, p + sig_blocks)
        if(length(mu_set) != nrow(Sigma)){
            print(i)
            stop("length(mu_set) != nrow(Sigma)")
        }
        X <- mvrnorm(n=n, mu=mu_set, Sigma=Sigma)
        
        # blocked_dgp_vars <- coefs[2]
        # sig_unblocked_vars <- coefs[3]
        mu <- as.numeric(X %*% beta)
        # print(beta)
        # Remove true blocked signal feature from x now that I've generated mu
        # (for now only have code written for one block)
        if(nblocks==1 & sig_blocks ==1){
            z <- X[, 1]
            X <- X[, 2:ncol(X)]
        } else{
            stop("!(nblocks==1 & sig_blocks ==1)")
        }
        # If SNR is null, use sigma_eps_sq
        if(!is.na(sigma_eps_sq)){
            sd <- sqrt(sigma_eps_sq)
        }else{
            sd <- sqrt(sum(mu^2) / (n * snr)) # taking snr = ||mu||^2 /(n * sigma^2)
        }
        
        y <- mu + sd * rnorm(n)
        ret_list[[i]] <- list(X=X, y=y)
    }
    
    return(ret_list)
}



random_simulate_func_ranking <- function(n, p, k_unblocked, beta_low_min,
    beta_low_max, beta_high, nblocks=1, sig_blocks=1, block_size, rho_low,
    rho_high, var, snr=NA, sigma_eps_sq=NA, Sigma, beta, nsim){
    if(is.na(snr) & is.na(sigma_eps_sq)){
        stop("Must specify one of snr or sigma_eps_sq")
    }

    # List we'll return: will be length nsim, and every element will be a
    # named list with elements X and y.
    ret_list <- list()
    stopifnot(nrow(Sigma) == p + sig_blocks & ncol(Sigma) == p + sig_blocks)
    if(!all(dim(Sigma) == c(p + sig_blocks, p + sig_blocks))){
        print("dim(Sigma):")
        print(dim(Sigma))
        print("p:")
        print(p)
        print("sig_blocks:")
        print(sig_blocks)
        stop("nrow(Sigma) != p + sig_blocks | ncol(Sigma) != p + sig_blocks")
    }
    for(i in 1:nsim){
        # print(Sigma)
        mu_set <- rep(0, p + sig_blocks)
        if(length(mu_set) != nrow(Sigma)){
            print(i)
            stop("length(mu_set) != nrow(Sigma)")
        }
        X <- mvrnorm(n=n, mu=mu_set, Sigma=Sigma)
        
        # blocked_dgp_vars <- coefs[2]
        # sig_unblocked_vars <- coefs[3]
        mu <- as.numeric(X %*% beta)
        # print(beta)
        # Remove true blocked signal feature from x now that I've generated mu
        # (for now only have code written for one block)
        if(nblocks==1 & sig_blocks==1){
            z <- X[, 1]
            X <- X[, 2:ncol(X)]
        } else{
            stop("!(nblocks==1 & sig_blocks==1)")
        }
        # If SNR is null, use sigma_eps_sq
        if(!is.na(sigma_eps_sq)){
            sd <- sqrt(sigma_eps_sq)
        }else{
            sd <- sqrt(sum(mu^2) / (n * snr)) # taking snr = ||mu||^2 /(n * sigma^2)
        }
        
        y <- mu + sd * rnorm(n)
        ret_list[[i]] <- list(X=X, y=y)
    }
    
    return(ret_list)
}



make_sparse_blocked_linear_model4_random <- function(n, p, k_unblocked, beta_low,
    beta_high, nblocks=1, sig_blocks=1, block_size, rho, var, snr=NA,
    sigma_eps_sq=NA) {
    # Same as make_sparse_blocked_linear_model, but sets blocked coefficients to
    # beta_high and unblocked coefficients to beta_low, and also removes the
    # true signal feature from x after generating the response y.
    if(is.na(snr) & is.na(sigma_eps_sq)){
        stop("Must specify one of snr or sigma_eps_sq")
    }

    Sigma <- make_covariance_matrix(p + sig_blocks, nblocks, block_size +
        sig_blocks, rho, var)
    coefs <- make_coefficients4(p + sig_blocks, k_unblocked, beta_low,
        beta_high, nblocks, sig_blocks, block_size + 1)
    beta <- coefs[[1]]
    
    my_model <- new_model(name = "sblm2_random", 
                label = sprintf("Linear model2 (random design) with correlated blocks (n = %s,
                    p = %s, k_unblocked = %s, nblocks = %s, sig_blocks = %s, rho = %s)",
                    n, p, k_unblocked, nblocks, sig_blocks, rho),
                params = list(n = n,
                    p = p, k_unblocked = k_unblocked, beta_low = beta_low,
                    beta_high = beta_high,
                    nblocks = nblocks, sig_blocks = sig_blocks,
                    block_size = block_size, rho = rho, var = var,
                    snr = snr,
                    sigma_eps_sq = sigma_eps_sq,
                    Sigma = Sigma,
                    beta = beta,
                    blocked_dgp_vars = coefs[[2]],
                    sig_unblocked_vars = coefs[[3]],
                    insig_blocked_vars = coefs[[4]]),
                    simulate = random_simulate_func
    )
    # print("finished generating model")
    return(my_model)
}




make_sparse_blocked_linear_model4_random_ranking <- function(n, p, k_unblocked,
    beta_low_min, beta_low_max, beta_high, nblocks=1, sig_blocks=1, block_size,
    rho_low, rho_high, var, snr=NA, sigma_eps_sq=NA) {
    # Same as make_sparse_blocked_linear_model_random, but ranges coefficients
    # of weak signal features from beta_low_min to beta_low_max in order to have
    # a definitive ranking of weak signal features, and similarly ranges
    # correlations of proxies with Z from rho_low to rho_high.
    if(is.na(snr) & is.na(sigma_eps_sq)){
        stop("Must specify one of snr or sigma_eps_sq")
    }

    Sigma <- make_covariance_matrix_ranking(p + sig_blocks, nblocks, block_size +
        sig_blocks, rho_low, rho_high, var)
    coefs <- make_coefficients4_ranking(p + sig_blocks, k_unblocked,
        beta_low_min, beta_low_max, beta_high, nblocks, sig_blocks,
        block_size + 1)
    beta <- coefs[[1]]
    
    my_model <- new_model(name = "sblm2_random_ranking", 
                label = sprintf("Linear model2 (random design) with correlated blocks and ranked features (n = %s,
                    p = %s, k_unblocked = %s, nblocks = %s, sig_blocks = %s, rho_low = %s, rho_high = %s)",
                    n, p, k_unblocked, nblocks, sig_blocks, rho_low, rho_high),
                params = list(n = n,
                    p = p, k_unblocked = k_unblocked,
                    beta_low_min = beta_low_min, beta_low_max = beta_low_max,
                    beta_high = beta_high,
                    nblocks = nblocks, sig_blocks = sig_blocks,
                    block_size = block_size, rho_low = rho_low,
                    rho_high=rho_high, var = var,
                    snr = snr,
                    sigma_eps_sq = sigma_eps_sq,
                    Sigma = Sigma,
                    beta = beta,
                    blocked_dgp_vars = coefs[[2]],
                    sig_unblocked_vars = coefs[[3]],
                    insig_blocked_vars = coefs[[4]]),
                    simulate = random_simulate_func_ranking
    )
    # print("finished generating model")
    return(my_model)
}




gen_mu_x_sd4_ranking <- function(n, p, k_unblocked, beta_low_min, beta_low_max,
    beta_high, nblocks=1, sig_blocks=1, block_size, rho_high, row_low, var,
    snr=NA, sigma_eps_sq=NA){
    if(is.na(snr) & is.na(sigma_eps_sq)){
        stop("Must specify one of snr or sigma_eps_sq")
    }
    Sigma <- make_covariance_matrix_ranking(p + sig_blocks, nblocks, block_size +
        sig_blocks, rho_low, rho_high, var)
    # print(Sigma)
    x <- mvrnorm(n=n, mu=rep(0, p + sig_blocks), Sigma=Sigma)
    coefs <- make_coefficients4_ranking(p + sig_blocks, k_unblocked,
        beta_low_min, beta_low_max, beta_high, nblocks, sig_blocks,
        block_size + 1)
    beta <- coefs[[1]]
    # blocked_dgp_vars <- coefs[2]
    # sig_unblocked_vars <- coefs[3]
    mu <- as.numeric(x %*% beta)
    # print(beta)
    # Remove true blocked signal feature from x now that I've generated mu
    # (for now only have code written for one block)
    if(nblocks==1 & sig_blocks ==1){
        z <- x[, 1]
        x <- x[, 2:ncol(x)]
    } else{
        stop("!(nblocks==1 & sig_blocks ==1)")
    }
    # If SNR is null, use sigma_eps_sq
    if(!is.na(sigma_eps_sq)){
        sd <- sqrt(sigma_eps_sq)
    }else{
        sd <- sqrt(sum(mu^2) / (n * snr)) # taking snr = ||mu||^2 /(n * sigma^2)
    }
    # Take out first row of Sigma that no longer needs to be there (again, 
    # assuming one block)
    if(nblocks==1 & sig_blocks ==1){
        Sigma <- Sigma[2:nrow(Sigma), 2:ncol(Sigma)]
    }
    return(list(mu, x, sd, Sigma, beta, z, coefs))
}
































gen_mu_x_sd4_weighted <- function(n, p, k_unblocked, beta_low, beta_high,
    nblocks=1, sig_blocks=1, block_size, n_strong_block_vars, rho_high, row_low,
    var, snr=NA, sigma_eps_sq=NA){

    # Generates a mu vector, design matrix x, standard deviation of noise
    # (calculated to ensure stated SNR), covariance matrix, coefficient vector,
    # z, coefficients

    if(is.na(snr) & is.na(sigma_eps_sq)){
        stop("Must specify one of snr or sigma_eps_sq")
    }
    Sigma <- make_covariance_matrix_weighted(p + sig_blocks, nblocks,
        block_size + sig_blocks, n_strong_block_vars + sig_blocks, rho_high,
        rho_low, var)

    coefs <- make_coefficients4_ranking2(p + sig_blocks, k_unblocked, beta_low,
        beta_high, nblocks, sig_blocks, block_size + 1)

    beta <- coefs$beta

    gen_mu_x_y_sd_res <- gen_mu_x_y_sd(n, p, beta, Sigma, sig_blocks,
        block_size, snr, sigma_eps_sq)

    x <- gen_mu_x_y_sd_res$X
    mu <- gen_mu_x_y_sd_res$mu
    blocked_dgp_vars <- gen_mu_x_y_sd_res$blocked_dgp_vars
    z <- gen_mu_x_y_sd_res$z
    sd <- gen_mu_x_y_sd_res$sd

    # Take out rows of Sigma corresponding to latent features

    Sigma <- Sigma[setdiff(1:ncol(Sigma), blocked_dgp_vars),
        setdiff(1:ncol(Sigma), blocked_dgp_vars)]

    # Confirm Sigma is square with dimension p
    if(nrow(Sigma) != ncol(x) | ncol(Sigma) != ncol(x)){
        # print("blocked_dgp_vars:")
        # print(blocked_dgp_vars)
        # print("dim(Sigma):")
        # print(dim(Sigma))
        # print("p:")
        # print(p)
        # print("dim(x):")
        # print(dim(x))
        stop("nrow(Sigma) != ncol(x) | ncol(Sigma) != ncol(x)")
    }

    # Confirm Sigma is symmetric
    if(any(Sigma != t(Sigma))){
        # print(Sigma)
        stop("any(Sigma != t(Sigma))")
    }

    # if(nblocks==1 & sig_blocks ==1){
    #     Sigma <- Sigma[2:nrow(Sigma), 2:ncol(Sigma)]
    # }
    return(list(mu=mu, x=x, sd=sd, Sigma=Sigma, beta=beta, z=z, coefs=coefs))
}

random_simulate_func_weighted <- function(n, n_clus, n_test, p, k_unblocked,
    est_clus_cutoff, beta_low, beta_high, nblocks=1, sig_blocks=1, block_size,
    n_strong_block_vars, rho_high, rho_low, var, snr=NA, sigma_eps_sq=NA, nsim){
    if(is.na(snr) & is.na(sigma_eps_sq)){
        stop("Must specify one of snr or sigma_eps_sq")
    }

    # List we'll return: will be length nsim, and every element will be a
    # named list with elements X, y, testX, testY, and testMu.
    ret_list <- list()

    for(i in 1:nsim){

        # Data for estimating clusters
        clus_data <- cssr::genClusteredDataWeighted(n=n_clus, p=p,
            k_unclustered=k_unblocked, cluster_size=block_size,
            n_strong_cluster_vars=n_strong_block_vars, n_clusters=nblocks,
            sig_clusters=sig_blocks, rho_high=rho_high, rho_low=rho_low,
            beta_latent=beta_high, beta_unclustered=beta_low, snr=snr,
            sigma_eps_sq=sigma_eps_sq)

        # Estimate clusters on clus_data
        est_clusters <- clustEsts(clus_data$X)

        gen_mu_x_y_sd_res <- cssr::genClusteredDataWeighted(n=n, p=p,
            k_unclustered=k_unblocked, cluster_size=block_size,
            n_strong_cluster_vars=n_strong_block_vars, n_clusters=nblocks,
            sig_clusters=sig_blocks, rho_high=rho_high, rho_low=rho_low,
            beta_latent=beta_high, beta_unclustered=beta_low, snr=snr,
            sigma_eps_sq=sigma_eps_sq)

        test_data <- cssr::genClusteredDataWeighted(n=n_test, p=p,
            k_unclustered=k_unblocked, cluster_size=block_size,
            n_strong_cluster_vars=n_strong_block_vars, n_clusters=nblocks,
            sig_clusters=sig_blocks, rho_high=rho_high, rho_low=rho_low,
            beta_latent=beta_high, beta_unclustered=beta_low, snr=snr,
            sigma_eps_sq=sigma_eps_sq)

        # All stability selection methods can use same output (other than
        # estimated vs. known clusters), so we will do this computation once
        # in the model generation and then use the results in different
        # ways depending on the method
        res_known <- getRes(gen_mu_x_y_sd_res$X, gen_mu_x_y_sd_res$y,
            clusters=1:block_size, nblocks=nblocks, sig_blocks=sig_blocks)

        # Confirm cluster in the results
        stopifnot(ncol(res_known$clus_sel_mat) == p - block_size + 1)

        res_est <- getRes(gen_mu_x_y_sd_res$X, gen_mu_x_y_sd_res$y,
            clusters=est_clusters, nblocks=nblocks, sig_blocks=sig_blocks)

        ret_list[[i]] <- list(est_clusters=est_clusters, X=gen_mu_x_y_sd_res$X,
            y=gen_mu_x_y_sd_res$y, res_known=res_known, res_est=res_est,
            testX=test_data$X, testY=test_data$y, testMu=test_data$mu)


    }
    
    return(ret_list)
}

sim_5_extra_sim_func <- function(n, n_clus, n_test, p, k_unblocked,
    est_clus_cutoff, beta_low, beta_high, nblocks=1, sig_blocks=1, block_size,
    rho_high, rho_low, var, snr=NA, sigma_eps_sq=NA, nsim){
    if(is.na(snr) & is.na(sigma_eps_sq)){
        stop("Must specify one of snr or sigma_eps_sq")
    }

    # List we'll return: will be length nsim, and every element will be a
    # named list with elements X, y, testX, testY, and testMu.
    ret_list <- list()

    for(i in 1:nsim){

        # Data for estimating clusters
        clus_data <- cssr::genClusteredDataWeightedRandom(n=n_clus, p=p,
            k_unclustered=k_unblocked, cluster_size=block_size,
            n_clusters=nblocks,
            sig_clusters=sig_blocks, rho_high=rho_high, rho_low=rho_low,
            beta_latent=beta_high, beta_unclustered=beta_low, snr=snr,
            sigma_eps_sq=sigma_eps_sq)

        # Estimate clusters on clus_data
        est_clusters <- clustEsts(clus_data$X)

        gen_mu_x_y_sd_res <- cssr::genClusteredDataWeightedRandom(n=n, p=p,
            k_unclustered=k_unblocked, cluster_size=block_size,
            n_clusters=nblocks,
            sig_clusters=sig_blocks, rho_high=rho_high, rho_low=rho_low,
            beta_latent=beta_high, beta_unclustered=beta_low, snr=snr,
            sigma_eps_sq=sigma_eps_sq)

        test_data <- cssr::genClusteredDataWeightedRandom(n=n_test, p=p,
            k_unclustered=k_unblocked, cluster_size=block_size,
            n_clusters=nblocks,
            sig_clusters=sig_blocks, rho_high=rho_high, rho_low=rho_low,
            beta_latent=beta_high, beta_unclustered=beta_low, snr=snr,
            sigma_eps_sq=sigma_eps_sq)

        # All stability selection methods can use same output (other than
        # estimated vs. known clusters), so we will do this computation once
        # in the model generation and then use the results in different
        # ways depending on the method
        res_known <- getRes(gen_mu_x_y_sd_res$X, gen_mu_x_y_sd_res$y,
            clusters=1:block_size, nblocks=nblocks, sig_blocks=sig_blocks)

        # Confirm cluster in the results
        stopifnot(ncol(res_known$clus_sel_mat) == p - block_size + 1)

        res_est <- getRes(gen_mu_x_y_sd_res$X, gen_mu_x_y_sd_res$y,
            clusters=est_clusters, nblocks=nblocks, sig_blocks=sig_blocks)

        ret_list[[i]] <- list(est_clusters=est_clusters, X=gen_mu_x_y_sd_res$X,
            y=gen_mu_x_y_sd_res$y, res_known=res_known, res_est=res_est,
            testX=test_data$X, testY=test_data$y, testMu=test_data$mu)
    }
    
    return(ret_list)
}



















make_covariance_matrix_weighted <- function(p, nblocks, block_size,
    n_strong_block_vars, rho_high, row_low, var) {
     # start with p x p identity matrix
    Sigma <- var*diag(p)

    ### select n.blocks blocks of features to be highly correlated 

    # create matrix with nblocks rows, each containing a vector of 
    # indices of highly correlated features
    block_feats <- matrix(seq(nblocks*block_size), nrow=nblocks, byrow=TRUE)

    # add covariances of highly correlated features to sigma
    for(i in 1:nblocks){
        for(j in 1:(block_size - 1)){
            for(k in (j+1):block_size){
                feat_1 <- block_feats[i, j]
                feat_2 <- block_feats[i, k]
                Sigma[feat_1, feat_2] <- rho_low
                Sigma[feat_2, feat_1] <- rho_low
            }
        }

        for(j in 1:(n_strong_block_vars - 1)){
            for(k in (j+1):n_strong_block_vars){
                feat_1 <- block_feats[i, j]
                feat_2 <- block_feats[i, k]
                Sigma[feat_1, feat_2] <- rho_high
                Sigma[feat_2, feat_1] <- rho_high
            }
        }
    }

    # Confirm Sigma is square with dimension p
    if(nrow(Sigma) != p | ncol(Sigma) != p){
        # print("dim(Sigma):")
        # print(dim(Sigma))
        # print("p:")
        # print(p)
        stop("nrow(Sigma) != p| ncol(Sigma) != p")
    }

    # Confirm Sigma is symmetric
    if(any(Sigma != t(Sigma))){
        # print(Sigma)
        stop("any(Sigma != t(Sigma))")
    }

    return(Sigma)
}

make_covariance_matrix_ranking <- function(p, nblocks, block_size,
    rho_low, row_high, var) {
     # start with p x p identity matrix
    Sigma <- var*diag(p)

    ### select n.blocks blocks of features to be highly correlated 

    # create matrix with nblocks rows, each containing a vector of 
    # indices of highly correlated features
    block_feats <- matrix(seq(nblocks*block_size), nrow=nblocks, byrow=TRUE)

    # Range of correlation coefficients
    rhos <- seq(rho_low, rho_high, length.out=block_size - 1)

    # add covariances of highly correlated features to sigma
   
    for(i in 1:nblocks){ 
        for(l in 1:block_size){
            j <- 1
            while(j < block_size - l + 1){
                if(block_size - l + 1 < j + 1){
                    stop("block_size - l + 1 < j + 1")
                }
                for(k in (j+1):(block_size - l + 1)){
                    feat_1 <- block_feats[i, j]
                    feat_2 <- block_feats[i, k]
                    # print("feat_1:")
                    # print(feat_1)
                    # print("feat_2:")
                    # print(feat_2)
                    Sigma[feat_1, feat_2] <- rhos[l]
                    Sigma[feat_2, feat_1] <- rhos[l]
                }
                j <- j + 1
            }
        }
    }



    stopifnot(nrow(Sigma) == p & ncol(Sigma) == p)

    stopifnot(all(Sigma == t(Sigma)))

    return(Sigma)
}








gen_mu_x_sd4_ranking2 <- function(n, p, k_unblocked, beta_low,
    beta_high, nblocks=1, sig_blocks=1, block_size, rho, var,
    snr=NA, sigma_eps_sq=NA){
    if(is.na(snr) & is.na(sigma_eps_sq)){
        stop("Must specify one of snr or sigma_eps_sq")
    }
    Sigma <- make_covariance_matrix(p + sig_blocks, nblocks, block_size +
        sig_blocks, rho, var)
    # print(Sigma)
    coefs <- make_coefficients4_ranking2(p + sig_blocks, k_unblocked,
        beta_low, beta_high, nblocks, sig_blocks, block_size + 1)
    beta <- coefs$beta


    gen_mu_x_y_sd_res <- gen_mu_x_y_sd(n, p, beta, Sigma, sig_blocks,
        block_size, snr, sigma_eps_sq)

    x <- gen_mu_x_y_sd_res$X
    mu <- gen_mu_x_y_sd_res$mu
    blocked_dgp_vars <- gen_mu_x_y_sd_res$blocked_dgp_vars
    z <- gen_mu_x_y_sd_res$z
    sd <- gen_mu_x_y_sd_res$sd


    # Take out first row of Sigma that no longer needs to be there (again, 
    # assuming one block)
    # Take out rows of Sigma corresponding to latent features

    Sigma <- Sigma[setdiff(1:ncol(Sigma), blocked_dgp_vars),
        setdiff(1:ncol(Sigma), blocked_dgp_vars)]

    # Confirm Sigma is square with dimension p
    if(nrow(Sigma) != ncol(x) | ncol(Sigma) != ncol(x)){
        # print("blocked_dgp_vars:")
        # print(blocked_dgp_vars)
        # print("dim(Sigma):")
        # print(dim(Sigma))
        # print("p:")
        # print(p)
        # print("dim(x):")
        # print(dim(x))
        stop("nrow(Sigma) != ncol(x) | ncol(Sigma) != ncol(x)")
    }

    # Confirm Sigma is symmetric
    if(any(Sigma != t(Sigma))){
        # print(Sigma)
        stop("any(Sigma != t(Sigma))")
    }

    return(list(mu=mu, x=x, sd=sd, Sigma=Sigma, beta=beta, z=z, coefs=coefs))
}

random_simulate_func_ranking2 <- function(n, n_clus, n_test, p, k_unblocked,
    est_clus_cutoff, beta_low, beta_high, nblocks=1, sig_blocks=1, block_size,
    rho, var, snr=NA, sigma_eps_sq=NA, nsim){
    if(is.na(snr) & is.na(sigma_eps_sq)){
        stop("Must specify one of snr or sigma_eps_sq")
    }

    # List we'll return: will be length nsim, and every element will be a
    # named list with elements, X, y, testX, testY, and testMu.

    ret_list <- list()

    for(i in 1:nsim){

        # Cluster estimation data
        clus_data <- cssr::genClusteredData(n=n_clus, p=p,
            k_unclustered=k_unblocked, cluster_size=block_size,
            n_clusters=nblocks, sig_clusters=sig_blocks, rho=rho,
            beta_latent=beta_high, beta_unclustered=beta_low, snr=snr,
            sigma_eps_sq=sigma_eps_sq)

        est_clusters <- clustEsts(clus_data$X)

        # Feature selection data
        gen_mu_x_y_sd_res <- cssr::genClusteredData(n=n, p=p,
            k_unclustered=k_unblocked, cluster_size=block_size,
            n_clusters=nblocks, sig_clusters=sig_blocks, rho=rho,
            beta_latent=beta_high, beta_unclustered=beta_low, snr=snr,
            sigma_eps_sq=sigma_eps_sq)

        # Training/test data
        test_data <- cssr::genClusteredData(n=n_test, p=p,
            k_unclustered=k_unblocked, cluster_size=block_size,
            n_clusters=nblocks, sig_clusters=sig_blocks, rho=rho,
            beta_latent=beta_high, beta_unclustered=beta_low, snr=snr,
            sigma_eps_sq=sigma_eps_sq)

        # All stability selection methods can use same output (other than
        # estimated vs. known clusters), so we will do this computation once
        # in the model generation and then use the results in different
        # ways depending on the method
        res_known <- getRes(gen_mu_x_y_sd_res$X, gen_mu_x_y_sd_res$y,
            clusters=1:block_size, nblocks=nblocks, sig_blocks=sig_blocks)

        # Confirm cluster in the results
        stopifnot(ncol(res_known$clus_sel_mat) == p - block_size + 1)

        res_est <- getRes(gen_mu_x_y_sd_res$X, gen_mu_x_y_sd_res$y,
            clusters=est_clusters, nblocks=nblocks, sig_blocks=sig_blocks)

        ret_list[[i]] <- list(est_clusters=est_clusters, X=gen_mu_x_y_sd_res$X,
            y=gen_mu_x_y_sd_res$y, res_known=res_known, res_est=res_est,
            testX=test_data$X, testY=test_data$y, testMu=test_data$mu)
    }
    
    return(ret_list)
}

getRes <- function(X, y, clusters, nblocks, sig_blocks){
    lambda <- cssr::getLassoLambda(X, y, lambda_choice="min")

    stopifnot(nblocks == 1)
    stopifnot(sig_blocks == 1)

    res <- cssr::css(X, y, lambda=lambda, clusters=clusters,
        num_cores=parallel::detectCores() - 1)

    return(res)
}

clustEsts <- function(X){
    # Estimate clusters on clus_data
    cor_mat <- stats::cor(X)
    dist <- stats::as.dist(1 - abs(cor_mat))
    h <- stats::hclust(dist)

    clus_assignments <-  stats::cutree(h, h=est_clus_cutoff)
    stopifnot(length(clus_assignments) == p)

    n_est_clusters <- length(unique(clus_assignments))

    est_clusters <- list()
    for(i in 1:n_est_clusters){
      est_clusters[[i]] <- which(clus_assignments == i)
    }

    return(est_clusters)
}

make_blocked_lin_mod4_ran_weight <- function(n, n_clus, n_test, p, k_unblocked, 
    beta_low, beta_high, nblocks=1, sig_blocks=1, block_size, est_clus_cutoff,
    n_strong_block_vars, rho_high, rho_low, var, snr=NA, sigma_eps_sq=NA) {
    # Same as make_sparse_blocked_linear_model4_random, but makes only
    # n_strong_block_vars have a high correlation with latent signal; remaining
    # block_size - n_strong_block_vars variables have low correlation.
    if(is.na(snr) & is.na(sigma_eps_sq)){
        stop("Must specify one of snr or sigma_eps_sq")
    }
    if(n_strong_block_vars < 0 | n_strong_block_vars > block_size){
        stop("n_strong_block_vars < 0 | n_strong_block_vars > 0 block_size")
    }

    if(rho_low >= rho_high){
        stop("rho_low >= rho_high")
    }

    # Make sure p is large enough
    stopifnot(p >= nblocks*block_size + k_unblocked)

    # Sigma <- make_covariance_matrix_weighted(p + sig_blocks, nblocks,
    #     block_size + sig_blocks, n_strong_block_vars + sig_blocks, rho_high,
    #     rho_low, var)

    # stopifnot(nrow(Sigma) == p + sig_blocks & ncol(Sigma) == p + sig_blocks)

    # coefs <- make_coefficients4_ranking2(p + sig_blocks, k_unblocked, beta_low,
    #     beta_high, nblocks, sig_blocks, block_size + 1)
    
    my_model <- new_model(name = "sblm2_random_weighted", 
                label = sprintf("Lin model (weight avg, corr blocks) (n= %s, p= %s, k_unblocked= %s, rho_high= %s)",
                    n, p, k_unblocked, rho_high),
                params = list(n = n, n_clus=n_clus, n_test = n_test,
                    p = p, k_unblocked = k_unblocked, beta_low = beta_low,
                    est_clus_cutoff=est_clus_cutoff, beta_high = beta_high,
                    nblocks = nblocks, sig_blocks = sig_blocks,
                    block_size = block_size, 
                    n_strong_block_vars = n_strong_block_vars,
                    rho_high = rho_high,
                    rho_low = rho_low, var = var,
                    snr = snr,
                    sigma_eps_sq = sigma_eps_sq
                    # Sigma = Sigma,
                    # beta = coefs$beta,
                    # blocked_dgp_vars = coefs$blocked_dgp_vars,
                    # sig_unblocked_vars = coefs$sig_unblocked_vars,
                    # insig_blocked_vars = coefs$insig_blocked_vars
                    )
                    , simulate = random_simulate_func_weighted
    )
    # print("finished generating model")
    return(my_model)
}

sim_5_extra <- function(n, n_clus, n_test, p, k_unblocked, 
    beta_low, beta_high, nblocks=1, sig_blocks=1, block_size, est_clus_cutoff,
    rho_high, rho_low, var, snr=NA, sigma_eps_sq=NA) {

    if(is.na(snr) & is.na(sigma_eps_sq)){
        stop("Must specify one of snr or sigma_eps_sq")
    }

    if(rho_low >= rho_high){
        stop("rho_low >= rho_high")
    }

    # Make sure p is large enough
    stopifnot(p >= nblocks*block_size + k_unblocked)
    
    my_model <- new_model(name = "sim_5_extra", 
                label = sprintf("Lin model (weight avg, corr blocks) (n= %s, p= %s, k_unblocked= %s, rho_high= %s)",
                    n, p, k_unblocked, rho_high),
                params = list(n = n, n_clus=n_clus, n_test = n_test,
                    p = p, k_unblocked = k_unblocked, beta_low = beta_low,
                    est_clus_cutoff=est_clus_cutoff, beta_high = beta_high,
                    nblocks = nblocks, sig_blocks = sig_blocks,
                    block_size = block_size, 
                    rho_high = rho_high,
                    rho_low = rho_low, var = var,
                    snr = snr,
                    sigma_eps_sq = sigma_eps_sq
                    )
                    , simulate = sim_5_extra_sim_func
    )
    # print("finished generating model")
    return(my_model)
}

make_sparse_blocked_linear_model4_random_ranking2 <- function(n, n_clus, n_test,
    p, k_unblocked, est_clus_cutoff, beta_low, beta_high, nblocks=1,
    sig_blocks=1, block_size, rho, var, snr=NA, sigma_eps_sq=NA) {
    # Same as make_sparse_blocked_linear_model_random, but ith coefficient
    # of weak signal features is beta_low/sqrt(i) in order to have
    # a definitive ranking of weak signal features.
    if(is.na(snr) & is.na(sigma_eps_sq)){
        stop("Must specify one of snr or sigma_eps_sq")
    }

    # Make sure p is large enough
    stopifnot(p >= nblocks*block_size + k_unblocked)

    # Sigma <- make_covariance_matrix(p + sig_blocks, nblocks, block_size +
    #     sig_blocks, rho, var)
    # coefs <- make_coefficients4_ranking2(p + sig_blocks, k_unblocked,
    #     beta_low, beta_high, nblocks, sig_blocks,
    #     block_size + 1)
    
    my_model <- new_model(name = "sblm2_random_ranking2", 
                label = sprintf("Linear model2 (random design) with correlated blocks and ranked2 features (n = %s,
                    p = %s, k_unblocked = %s, nblocks = %s, sig_blocks = %s, rho = %s)",
                    n, p, k_unblocked, nblocks, sig_blocks, rho),
                params = list(n = n, n_clus=n_clus, n_test = n_test,
                    p = p, k_unblocked = k_unblocked,
                    est_clus_cutoff=est_clus_cutoff, beta_low = beta_low, 
                    beta_high = beta_high,
                    nblocks = nblocks, sig_blocks = sig_blocks,
                    block_size = block_size, rho= rho, var = var,
                    snr = snr,
                    sigma_eps_sq = sigma_eps_sq
                    # blocked_dgp_vars = coefs$blocked_dgp_vars,
                    # sig_unblocked_vars = coefs$sig_unblocked_vars,
                    # insig_blocked_vars = coefs$insig_blocked_vars
                    )
                    , simulate = random_simulate_func_ranking2
    )
    # print("finished generating model")
    return(my_model)
}

plant_model <- function(n_selec, n_train, n_test, cor_cutoff, response_name,
    max_model_size) {
    
    plant_mod <- new_model(name = "plant_mod"
        , label = "GWAS Simulation Study Model"
        , params = list(n_selec=n_selec, n_train=n_train, n_test=n_test,
            cor_cutoff=cor_cutoff, response_name=response_name,
            max_model_size=max_model_size)
        , simulate = plant_sim_func
    )
    return(plant_mod)
}

plant_sim_func <- function(n_selec, n_train, n_test, cor_cutoff,
    response_name, max_model_size, nsim){
    n <- n_selec + n_train + n_test
    stopifnot(n == nrow(snps))
    stopifnot(n == nrow(data_plant))
    # List we'll return: will be length nsim, and every element will be a
    # named list with selec_inds, train_inds, and test_inds

    ret_list <- list()

    for(i in 1:nsim){

        selec_inds <- sample.int(n, size=n_selec, replace=FALSE)
        train_inds <- sample(setdiff(1:n, selec_inds), size=n_train,
            replace=FALSE)
        test_inds <- setdiff(1:n, c(selec_inds, train_inds))

        stopifnot(length(test_inds) == n_test)

        # Get clusters from non-test set data
        est_clusters <- clustEstsPlant(snps_mat=snps[c(selec_inds, train_inds),
            ], cor_cutoff=cor_cutoff)

        # Get cluster stability selection proportions on selection data
        selec_data <- getPlantSelecData(selec_inds, response_name)
        X_selec <- selec_data$X
        y_selec <- selec_data$y
        res_est <- getResPlant(X_selec, y_selec, clusters=est_clusters)

        ret_list[[i]] <- list(selec_inds=selec_inds, train_inds=train_inds,
            test_inds=test_inds, est_clusters=est_clusters, res_est=res_est)
    }
    
    return(ret_list)
}

getResPlant <- function(X, y, clusters){
    lambda <- cssr::getLassoLambda(X, y, lambda_choice="min")

    res <- cssr::css(X, y, lambda=lambda, clusters=clusters,
        num_cores=parallel::detectCores() - 1)

    return(res)
}

clustEstsPlant <- function(snps_mat, cor_cutoff){
    cor_mat <- cor(snps_mat)

    zero_sds <- apply(snps_mat, 2, function(x){length(unique(x)) == 1})

    # Should have already confirmed there are no SNPs with zero standard
    # deviation
    stopifnot(sum(zero_sds) == 0)

    # zero_sd_inds <- which(zero_sds)

    # for(i in zero_sd_inds){
    #     cor_mat[, i] <- 0
    #     cor_mat[i, ] <- 0
    #     cor_mat[i, i] <- 1
    # }

    dist <- stats::as.dist(1 - abs(cor_mat))
    h <- stats::hclust(dist)

    clus_assignments <-  stats::cutree(h, h=cor_cutoff)
    stopifnot(length(clus_assignments) == n_snps)

    n_est_clusters <- length(unique(clus_assignments))

    est_clusters <- list()
    for(i in 1:n_est_clusters){
      est_clusters[[i]] <- which(clus_assignments == i)
    }

    return(est_clusters)
}