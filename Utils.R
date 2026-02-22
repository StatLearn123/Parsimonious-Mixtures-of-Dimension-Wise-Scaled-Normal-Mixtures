# library(mclust) # must be installed and exported in foreach
# library(mnormt) # must be installed
# library(expint) # must be installed
# library(SphericalCubature) # must be installed

#### Scale-correlation decomposition ####

scale_correlation_decomp <- function(Sigma) {
  # Extract standard deviations (square root of diagonal elements)
  D <- sqrt(diag(Sigma))

  # Construct the diagonal matrix
  D_matrix <- diag(D)

  # Compute the correlation matrix
  R <- solve(D_matrix) %*% Sigma %*% solve(D_matrix)

  return(list(T = D_matrix, P = R))
}

#### Number of parameters ####

npar.pars <- function(k, d, corr = "Free", scale = "Free", tailedness = "Free") {
  # cor: Free, E.comp, Identity
  # scale: Free, E.comp, E.dim, E.comp.dim
  # tailedness: Free, E.comp, E.dim, E.comp.dim, Normal

  npar <- (k - 1) + k * d

  if (corr == "Free") {
    npar <- npar + k * d * (d - 1) / 2
  }
  if (corr == "E.comp") {
    npar <- npar + d * (d - 1) / 2
  }

  if (scale == "Free") {
    npar <- npar + k * d
  }
  if (scale == "E.comp") {
    npar <- npar + d
  }
  if (scale == "E.dim") {
    npar <- npar + k
  }
  if (scale == "E.comp.dim") {
    npar <- npar + 1
  }

  if (tailedness == "Free") {
    npar <- npar + k * d
  }
  if (tailedness == "E.comp") {
    npar <- npar + d
  }
  if (tailedness == "E.dim") {
    npar <- npar + k
  }
  if (tailedness == "E.comp.dim") {
    npar <- npar + 1
  }

  return(npar)
}

#### Fit, via EM-based algorithms, of the DSNM mixture (no parallel) ####

DSMN.clust.new <- function(
    X, # data matrix
    k, # number of mixture components
    corr = "Free",
    scale = "Free",
    tailedness = "Free",
    rel.tol = 0.001, # one of the 2 stopping rules (see iter.max for the other)
    iter.max = 1000, # maximum number of iterations
    method = "L-BFGS-B",
    method2 = "Nelder-Mead",
    Asympt = 20) {
  if (is.vector(X)) {
    X <- matrix(X, ncol = 1)
  }
  if (is.data.frame(X)) {
    X <- as.matrix(X)
  }
  if (!is.matrix(X)) {
    X <- matrix(X, nrow = length(X), ncol = 1, dimnames = list(names(X), deparse(substitute(X))))
  }
  if (!is.numeric(X)) {
    stop("numeric matrix/vector expected for X")
  }
  if (any(is.na(X))) {
    stop("No NAs allowed")
  }
  if (is.null(k)) {
    stop("The number of mixture components k is NULL")
  }

  ptm <- proc.time()
  name <- as.character(c(corr, scale, tailedness))

  n <- nrow(X) # sample size
  d <- ncol(X) # number of variables

  # --------------------- #
  # Definition of objects #
  # --------------------- #

  # E-step quantities

  z <- densX <- dens <- array(NA, c(n, k), dimnames = list(1:n, paste("comp.", 1:k, sep = "")))
  w <- array(2, c(n, d, k), dimnames = list(1:n, paste("dim.", 1:d, sep = ""), paste("comp.", 1:k, sep = "")))
  delta <- array(NA, c(n, d, k), dimnames = list(1:n, paste("dim.", 1:d, sep = ""), paste("comp.", 1:k, sep = "")))

  # M-step quantities; sigma are the square root of the diagonal values of Sigma

  mu <- sigma <- array(NA, c(d, k), dimnames = list(paste("dim.", 1:d, sep = ""), paste("comp.", 1:k, sep = "")))
  Sigma <- P <- array(NA, c(d, d, k), dimnames = list(paste("dim.", 1:d, sep = ""), paste("dim.", 1:d, sep = ""), paste("comp.", 1:k, sep = "")))
  mu.rep <- array(NA, c(n, d, k), dimnames = list(paste("unit.", 1:n, sep = ""), paste("dim.", 1:d, sep = ""), paste("comp.", 1:k, sep = "")))
  sigma.rep <- array(NA, c(n, d, k), dimnames = list(paste("unit.", 1:n, sep = ""), paste("dim.", 1:d, sep = ""), paste("comp.", 1:k, sep = "")))
  theta <- array(NA, c(d, k), dimnames = list(paste("dim.", 1:d, sep = ""), paste("comp.", 1:k, sep = "")))
  tau_values <- NA

  # ------------------- #
  # Initialization of z #
  # ------------------- #

  temp <- mclust::Mclust(data = X, G = k, modelNames = "VVV", verbose = F)
  z <- temp$z

  if (corr != "Identity" | corr != scale) {
    for (g in 1:k) {
      temp2 <- cov.wt(x = X, wt = z[, g], cor = TRUE, center = TRUE, method = "ML")
      Sigma[, , g] <- temp2$cov
      prior <- colMeans(z)

      if (d > 1) {
        sigma[, g] <- sqrt(diag(Sigma[, , g]))
      } else {
        sigma[, g] <- sqrt(Sigma[, , g])
      }
    }

    # E.comp, E.dim, and E.comp.dim

    if (scale == "E.comp") {
      temp.sigma <- sqrt(rowMeans(sigma^2))
      for (g in 1:k) {
        sigma[, g] <- temp.sigma
      }
    }
    if (scale == "E.dim") {
      temp.sigma <- sqrt(colMeans(sigma^2))
      for (j in 1:d) {
        sigma[j, ] <- temp.sigma
      }
    }
    if (scale == "E.comp.dim") {
      temp.sigma <- sqrt(mean(sigma^2))
      sigma <- array(rep(temp.sigma, d * k), c(d, k), dimnames = list(paste("dim.", 1:d, sep = ""), paste("comp.", 1:k, sep = "")))
    }

    for (g in 1:k) {
      sigma.rep[, , g] <- matrix(rep(sigma[, g], n), n, d, byrow = TRUE)
    }
  }

  if (corr == "Identity") {
    for (g in 1:k) {
      P[, , g] <- diag(d)
    }
  } else if (corr == "E.comp") {
    for (g in 1:k) {
      P[, , g] <- cov2cor(apply(Sigma, 1:2, weighted.mean, w = prior))
    }
  } else if (corr == "Free") {
    for (g in 1:k) {
      P[, , g] <- cov2cor(Sigma[, , g])
    }
  }

  # Preliminary definition of the convergence criteria

  loglik <- -Inf
  check <- mark <- 0
  iteration <- 1

  while (check < 1) {
    cat("*")

    # ++++++ #
    # M-Step #
    # ++++++ #

    # .. #
    # pi #
    # .. #

    sizes <- colSums(z)
    prior <- colMeans(z)

    # .. #
    # mu #
    # .. #

    for (g in 1:k) {
      tempM1 <- solve(Reduce(`+`, lapply(1:n, function(i) z[i, g] * diag(sqrt(w[i, , g])) %*% solve(Sigma[, , g]) %*% diag(sqrt(w[i, , g])))))
      tempM2 <- Reduce(`+`, lapply(1:n, function(i) z[i, g] * diag(sqrt(w[i, , g])) %*% solve(Sigma[, , g]) %*% diag(sqrt(w[i, , g])) %*% X[i, ]))
      mu[, g] <- c(tempM1 %*% tempM2)

      mu.rep[, , g] <- matrix(rep(mu[, g], n), n, d, byrow = TRUE)
    }

    # .......... #
    # Sigma & Co #
    # .......... #

    if (corr == "Identity") {
      tau_values <- update_tau_I(x = X, z = z, w = w, mu = mu, model_type = scale)
    }

    if (corr == "E.comp") {
      res <- update_tau_P(
        x = X, z = z, w = w, mu = mu, model_type = scale, Sigma = Sigma, prior = prior, sigma = sigma, P = P,
        method = method, method2 = method2, maxit = 1000, trace = 0
      )

      tau_values <- res$tau
      P <- res$P
    }

    if (corr == "Free") {
      res <- update_tau_Pg(
        x = X, z = z, w = w, mu = mu, model_type = scale, Sigma = Sigma, prior = prior, sigma = sigma, P = P,
        method = method, method2 = method2, maxit = 1000, trace = 0
      )

      tau_values <- res$tau
      P <- res$P
    }

    for (g in 1:k) {
      Sigma[, , g] <- tau_values[, , g] %*% P[, , g] %*% tau_values[, , g]
      sigma[, g] <- diag(tau_values[, , g])
      sigma.rep[, , g] <- matrix(rep(sigma[, g], n), n, d, byrow = TRUE)
    }

    # ..... #
    # theta #
    # ..... #

    theta <- update_theta(z = z, w = w, n = n, d = d, k = k, model_type = tailedness, Asympt = Asympt)

    # +++++++++++++++++++ #
    # Density Computation #
    # +++++++++++++++++++ #

    for (g in 1:k) {
      densX[, g] <- dDSSEN(X, mu = mu[, g], Sigma = Sigma[, , g], theta = theta[, g])
      dens[, g] <- prior[g] * densX[, g]
    }
    dens.mixt <- rowSums(dens)

    # ------------------------------------- #
    # Global - Observed-data log-likelihood #
    # ------------------------------------- #

    llvalues <- sum(log(dens.mixt))
    loglik <- c(loglik, llvalues)

    iteration <- iteration + 1

    if (loglik[iteration] < loglik[iteration - 1]) {
      mark <- 1
    }

    if (iteration == iter.max | (loglik[iteration] - loglik[iteration - 1]) < rel.tol) {
      check <- 1
    }

    # ++++++ #
    # E-step #
    # ++++++ #

    # posterior probabilities of group membership

    z <- dens / matrix(rep(dens.mixt, k), ncol = k)

    for (g in 1:k) {
      delta[, , g] <- (X - mu.rep[, , g])^2 / sigma.rep[, , g]^2
      for (i in 1:n) {
        for (h in 1:d) {
          num <- expint::gammainc(a = (1 / 2 + 2), x = (delta[i, h, g] / 2 + theta[h, g]))
          den <- (delta[i, h, g] / 2 + theta[h, g]) * expint::gammainc(a = (1 / 2 + 1), x = (delta[i, h, g] / 2 + theta[h, g]))
          w[i, h, g] <- num / den

          if (is.nan(w[i, h, g])) w[i, h, g] <- 1
          if (w[i, h, g] == Inf) w[i, h, g] <- 50
          if (w[i, h, g] < 1) w[i, h, g] <- 1.000001
        }
      }
    }
  }

  loglik.final <- loglik[iteration]

  # The EM-algorithm is finished #

  # --------------------- #
  # Classification Matrix #
  # --------------------- #

  clust.vec <- mclust::map(z)
  clust.mat <- mclust::unmap(z)

  # -------------------- #
  # Number of parameters #
  # -------------------- #

  npar <- npar.pars(k = k, d = d, corr = corr, scale = scale, tailedness = tailedness)

  # -------------------- #
  # Information criteria #
  # -------------------- #

  BIC <- -2 * loglik.final + npar * log(n) # to be minimized

  z.const <- (z < 10^(-322)) * 10^(-322) + (z > 10^(-322)) * z
  hard.z <- (matrix(rep(apply(z, 1, max), k), n, k, byrow = F) == z) * 1

  EN <- -sum(hard.z * log(z.const))
  ICL <- BIC + 2 * EN

  ptm2 <- proc.time() - ptm
  time <- ptm2[3]

  result <- list(
    name           = name,
    time           = time,
    dens.comp      = dens,
    dens.mixt      = dens.mixt,
    z              = z,
    clust.mat      = clust.mat,
    clust.vec      = clust.vec,
    iter.stop      = iteration,
    loglik         = loglik.final,
    mark           = mark,
    BIC            = BIC,
    ICL            = ICL,
    sizes          = sizes,
    prior          = prior,
    mu             = mu,
    P              = P,
    tau_values     = tau_values,
    sigma          = sigma,
    Sigma          = Sigma,
    theta          = theta
  )
  class(result) <- "DSNM.mixt"

  return(result)
}

#### Function to update tau ####

update_tau_I <- function(x, z, w, mu, model_type) {
  n <- nrow(x) # Number of observations
  d <- ncol(x) # Number of dimensions
  k <- dim(z)[2] # Number of clusters

  tau <- array(0, dim = c(d, d, k)) # Placeholder for tau values
  temp <- array(NA, dim = c(n, d, k))

  for (i in 1:n) {
    for (j in 1:d) {
      for (g in 1:k) {
        temp[i, j, g] <- z[i, g] * w[i, j, g] * (x[i, j] - mu[j, g])^2
      }
    }
  }

  if (model_type == "E.comp.dim") {
    tau_value <- sqrt(sum(temp) / (n * d))

    for (g in 1:k) {
      tau[, , g] <- diag(tau_value, d)
    }
  } else if (model_type == "E.dim") {
    for (g in 1:k) {
      tau_value <- sqrt(sum(temp[, , g]) / (sum(z[, g]) * d))
      tau[, , g] <- diag(tau_value, d)
    }
  } else if (model_type == "E.comp") {
    for (j in 1:d) {
      tau_value <- sqrt(sum(temp[, j, ]) / n)
      tau[j, j, ] <- tau_value
    }
  } else if (model_type == "Free") {
    for (j in 1:d) {
      for (g in 1:k) {
        tau_value <- sqrt(sum(temp[, j, g]) / sum(z[, g]))
        tau[j, j, g] <- tau_value
      }
    }
  }

  return(tau)
}

#### Functions to update tau and P - Part I ####

obj_fun_P_E.comp.dim <- function(par, C, Z, tau_sq, n) {
  P <- l2P(par)

  P <- (P + t(P)) / 2

  term1 <- -(n / 2) * log(det(P))
  term2 <- -(1 / (2 * tau_sq^2)) * sum(diag(solve(P) %*% C %*% Z %*% t(C)))

  return(term1 + term2)
}

obj_fun_P_E.dim <- function(par, C, A, n) {
  P <- l2P(par)

  P <- (P + t(P)) / 2

  term1 <- -(n / 2) * log(det(P))
  term2 <- -(1 / 2) * sum(diag(solve(P) %*% C %*% A %*% t(C)))

  return(term1 + term2)
}

obj_fun_P_Free <- function(par, B, Z, n) {
  P <- l2P(par)

  P <- (P + t(P)) / 2

  term1 <- -(n / 2) * log(det(P))
  term2 <- -(1 / 2) * sum(diag(solve(P) %*% B %*% Z %*% t(B)))

  return(term1 + term2)
}

obj_fun_P_Free_tau <- function(par, z, S, d) {
  T1 <- exp(par)
  T <- diag(T1)

  term1 <- -sum(z) * log(det(T))
  term2 <- -(1 / 2) * t(rep(1, d)) %*% solve(T) %*% S %*% solve(T) %*% (rep(1, d))

  return(term1 + term2)
}

update_tau_P <- function(x, z, w, mu, model_type, Sigma, prior, sigma, P, method = "L-BFGS-B", method2 = "Nelder-Mead", maxit = 2000, trace = 0, lower = -4, upper = 4) {
  n <- nrow(x) # Number of observations
  d <- ncol(x) # Number of dimensions
  k <- dim(z)[2] # Number of clusters

  tau <- P_fin <- array(NA, dim = c(d, d, k)) # Placeholder for tau values and P
  Z <- diag(as.vector(z))
  C <- matrix(NA, d, n * k)

  P_init <- P2l(P[, , 1])

  if (model_type == "E.comp.dim") {
    cnt <- 0
    for (g in 1:k) {
      for (i in 1:n) {
        cnt <- cnt + 1
        C[, cnt] <- diag(sqrt(w[i, , g])) %*% (X[i, ] - mu[, g])
      }
    }

    if (method == "L-BFGS-B") {
      opt_res <- optim(
        par = P_init, fn = obj_fun_P_E.comp.dim, method = method, control = list(fnscale = -1, maxit = maxit, trace = trace),
        C = C, Z = Z, tau_sq = mean(sigma), n = n, lower = lower, upper = upper
      )
    } else {
      opt_res <- optim(
        par = P_init, fn = obj_fun_P_E.comp.dim, method = method, control = list(fnscale = -1, maxit = maxit, trace = trace),
        C = C, Z = Z, tau_sq = mean(sigma), n = n
      )
    }

    P_opt <- l2P(opt_res$par)
    P_est <- (P_opt + t(P_opt)) / 2

    tau_value <- sqrt(sum(diag(solve(P_est) %*% C %*% Z %*% t(C))) / (n * d))
    for (g in 1:k) {
      tau[, , g] <- diag(tau_value, d)
      P_fin[, , g] <- P_est
    }
  } else if (model_type == "E.dim") {
    A <- kronecker(diag(sigma[1, ]^(-2)), diag(n)) %*% Z

    cnt <- 0
    for (g in 1:k) {
      for (i in 1:n) {
        cnt <- cnt + 1
        C[, cnt] <- diag(sqrt(w[i, , g])) %*% (X[i, ] - mu[, g])
      }
    }

    if (method == "L-BFGS-B") {
      opt_res <- optim(
        par = P_init, fn = obj_fun_P_E.dim, method = method, control = list(fnscale = -1, maxit = maxit, trace = trace),
        C = C, A = A, n = n, lower = lower, upper = upper
      )
    } else {
      opt_res <- optim(
        par = P_init, fn = obj_fun_P_E.dim, method = method, control = list(fnscale = -1, maxit = maxit, trace = trace),
        C = C, A = A, n = n
      )
    }

    P_opt <- l2P(opt_res$par)
    P_est <- (P_opt + t(P_opt)) / 2

    for (g in 1:k) {
      C_g <- C[, ((g - 1) * n + 1):(g * n)]
      Z_g <- diag(z[, g])

      tau_value <- sqrt(sum(diag(solve(P_est) %*% C_g %*% Z_g %*% t(C_g))) / (sum(z[, g]) * d))
      tau[, , g] <- diag(tau_value, d)
      P_fin[, , g] <- P_est
    }
  } else if (model_type == "E.comp") {
    c <- array(NA, dim = c(n, d, k))

    for (g in 1:k) {
      for (i in 1:n) {
        c[i, , g] <- diag(sqrt(w[i, , g])) %*% (X[i, ] - mu[, g])
      }
    }

    temp <- array(NA, dim = c(d, d, n, k))

    for (i in 1:n) {
      for (g in 1:k) {
        temp[, , i, g] <- z[i, g] * c[i, , g] %*% t(c[i, , g])
      }
    }

    Sig <- rowSums(temp, dims = 2) / n
    Sig.dec <- scale_correlation_decomp(Sig)

    for (g in 1:k) {
      tau[, , g] <- Sig.dec$T
      P_fin[, , g] <- Sig.dec$P
    }
  } else if (model_type == "Free") {
    ### P ###

    cnt <- 0
    for (g in 1:k) {
      for (i in 1:n) {
        cnt <- cnt + 1
        C[, cnt] <- diag(sqrt(w[i, , g])) %*% (X[i, ] - mu[, g])
      }
    }

    B_list <- list()
    for (g in 1:k) {
      C_g <- C[, ((g - 1) * n + 1):(g * n)]

      B_list[[g]] <- solve(diag(sigma[, g])) %*% C_g
    }

    B <- do.call(cbind, B_list)

    if (method == "L-BFGS-B") {
      opt_res <- optim(
        par = P_init, fn = obj_fun_P_Free, method = method, control = list(fnscale = -1, maxit = maxit, trace = trace),
        B = B, Z = Z, n = n, lower = lower, upper = upper
      )
    } else {
      opt_res <- optim(
        par = P_init, fn = obj_fun_P_Free, method = method, control = list(fnscale = -1, maxit = maxit, trace = trace),
        B = B, Z = Z, n = n
      )
    }

    P_opt <- l2P(opt_res$par)
    P_est <- (P_opt + t(P_opt)) / 2

    ### T ###

    S_g <- array(NA, dim = c(d, d, n))
    c <- array(NA, dim = c(n, d, k))

    for (g in 1:k) {
      for (i in 1:n) {
        c[i, , g] <- diag(sqrt(w[i, , g])) %*% (X[i, ] - mu[, g])
        S_g[, , i] <- z[i, g] * diag(c[i, , g]) %*% solve(P_est) %*% diag(c[i, , g])
      }

      T_init <- log(sigma[, g])

      opt_res2 <- optim(
        par = T_init, fn = obj_fun_P_Free_tau, method = method2, control = list(fnscale = -1, maxit = maxit, trace = trace),
        z = z[, g], S = rowSums(S_g, dims = 2), d = d
      )

      # Revert the transformation
      T_opt <- exp(opt_res2$par)

      # Get the optimized parameters and transform them back into the matrix
      tau[, , g] <- diag(T_opt, d)
      P_fin[, , g] <- P_est
    }
  }

  return(list(P = P_fin, tau = tau))
}

#### Functions to update tau and P - Part II ####

obj_fun_Pg_E.comp.dim <- function(par, C, Z, tau_sq, d, k, z) {
  P2 <- matrix(par, d * (d - 1) / 2, k)
  P <- array(NA, dim = c(d, d, k))

  for (g in 1:k) {
    P[, , g] <- l2P(P2[, g])

    P[, , g] <- (P[, , g] + t(P[, , g])) / 2
  }

  term1 <- term2 <- numeric(k)

  for (g in 1:k) {
    term1[g] <- sum(z[, g]) * log(det(P[, , g]))

    term2[g] <- sum(diag(solve(P[, , g]) %*% C[, , g] %*% Z[, , g] %*% t(C[, , g])))
  }

  term3 <- -(1 / 2) * sum(term1) - (1 / (2 * tau_sq^2)) * sum(term2)

  return(term3)
}

obj_fun_Pg_E.dim <- function(par, C, Z, tau_sq, d, k, z) {
  P2 <- matrix(par, d * (d - 1) / 2, k)
  P <- array(NA, dim = c(d, d, k))

  for (g in 1:k) {
    P[, , g] <- l2P(P2[, g])

    P[, , g] <- (P[, , g] + t(P[, , g])) / 2
  }

  term1 <- term2 <- numeric(k)

  for (g in 1:k) {
    term1[g] <- sum(z[, g]) * log(det(P[, , g]))

    term2[g] <- sum(diag(solve(P[, , g]) %*% C[, , g] %*% Z[, , g] %*% t(C[, , g]))) / (tau_sq[g]^2)
  }

  term3 <- -(1 / 2) * sum(term1) - (1 / 2) * sum(term2)

  return(term3)
}

obj_fun_Pg_E.comp <- function(par, M, Z, d, k, z) {
  P2 <- matrix(par, d * (d - 1) / 2, k)
  P <- array(NA, dim = c(d, d, k))

  for (g in 1:k) {
    P[, , g] <- l2P(P2[, g])

    P[, , g] <- (P[, , g] + t(P[, , g])) / 2
  }

  term1 <- term2 <- numeric(k)

  for (g in 1:k) {
    term1[g] <- sum(z[, g]) * log(det(P[, , g]))

    term2[g] <- sum(diag(solve(P[, , g]) %*% M[, , g] %*% Z[, , g] %*% t(M[, , g])))
  }

  term3 <- -(1 / 2) * sum(term1) - (1 / 2) * sum(term2)

  return(term3)
}

obj_fun_Pg_E.comp_tau <- function(par, n, S, d) {
  T1 <- exp(par)
  T <- diag(T1)

  term1 <- -n * log(det(T))
  term2 <- -(1 / 2) * t(rep(1, d)) %*% solve(T) %*% S %*% solve(T) %*% (rep(1, d))

  return(term1 + term2)
}

update_tau_Pg <- function(x, z, w, mu, model_type, Sigma, prior, sigma, P, method = "L-BFGS-B", method2 = "Nelder-Mead", maxit = 2000, trace = 0, lower = -4, upper = 4) {
  n <- nrow(x) # Number of observations
  d <- ncol(x) # Number of dimensions
  k <- dim(z)[2] # Number of clusters

  tau <- P_fin <- array(NA, dim = c(d, d, k)) # Placeholder for tau values and P

  Z_g <- array(NA, dim = c(n, n, k))
  C_g <- array(NA, dim = c(d, n, k))
  P_init <- matrix(NA, nrow = d * (d - 1) / 2, ncol = k)

  if (model_type == "E.comp.dim") {
    for (g in 1:k) {
      Z_g[, , g] <- diag(z[, g])

      for (i in 1:n) {
        C_g[, i, g] <- diag(sqrt(w[i, , g])) %*% (X[i, ] - mu[, g])
      }

      P_init[, g] <- P2l(P[, , g])
    }

    P_init <- as.vector(P_init)

    if (method == "L-BFGS-B") {
      opt_res <- optim(
        par = P_init, fn = obj_fun_Pg_E.comp.dim, method = method, control = list(fnscale = -1, maxit = maxit, trace = trace),
        C = C_g, Z = Z_g, tau_sq = mean(sigma), d = d, k = k, z = z, lower = lower, upper = upper
      )
    } else {
      opt_res <- optim(
        par = P_init, fn = obj_fun_Pg_E.comp.dim, method = method, control = list(fnscale = -1, maxit = maxit, trace = trace),
        C = C_g, Z = Z_g, tau_sq = mean(sigma), d = d, k = k, z = z
      )
    }

    P_opt2 <- matrix(opt_res$par, d * (d - 1) / 2, k)
    P_est <- array(NA, dim = c(d, d, k))

    for (g in 1:k) {
      P_est[, , g] <- l2P(P_opt2[, g])
      P_est[, , g] <- (P_est[, , g] + t(P_est[, , g])) / 2
    }

    tau_value0 <- numeric(k)

    for (g in 1:k) {
      tau_value0[g] <- sum(diag(solve(P_est[, , g]) %*% C_g[, , g] %*% Z_g[, , g] %*% t(C_g[, , g])))
    }

    tau_value <- sqrt(sum(tau_value0) / (n * d))

    for (g in 1:k) {
      tau[, , g] <- diag(tau_value, d)
      P_fin[, , g] <- P_est[, , g]
    }
  } else if (model_type == "E.dim") {
    for (g in 1:k) {
      Z_g[, , g] <- diag(z[, g])

      for (i in 1:n) {
        C_g[, i, g] <- diag(sqrt(w[i, , g])) %*% (X[i, ] - mu[, g])
      }

      P_init[, g] <- P2l(P[, , g])
    }

    P_init <- as.vector(P_init)

    if (method == "L-BFGS-B") {
      opt_res <- optim(
        par = P_init, fn = obj_fun_Pg_E.dim, method = method, control = list(fnscale = -1, maxit = maxit, trace = trace),
        C = C_g, Z = Z_g, tau_sq = sigma[1, ], d = d, k = k, z = z, lower = lower, upper = upper
      )
    } else {
      opt_res <- optim(
        par = P_init, fn = obj_fun_Pg_E.dim, method = method, control = list(fnscale = -1, maxit = maxit, trace = trace),
        C = C_g, Z = Z_g, tau_sq = sigma[1, ], d = d, k = k, z = z
      )
    }

    P_opt2 <- matrix(opt_res$par, d * (d - 1) / 2, k)
    P_est <- array(NA, dim = c(d, d, k))

    for (g in 1:k) {
      P_est[, , g] <- l2P(P_opt2[, g])
      P_est[, , g] <- (P_est[, , g] + t(P_est[, , g])) / 2
    }

    temp <- matrix(NA, n, k)
    c <- array(NA, dim = c(n, d, k))

    for (g in 1:k) {
      for (i in 1:n) {
        c[i, , g] <- diag(sqrt(w[i, , g])) %*% (X[i, ] - mu[, g])
      }
    }

    for (i in 1:n) {
      for (g in 1:k) {
        temp[i, g] <- z[i, g] * t(c[i, , g]) %*% solve(P_est[, , g]) %*% c[i, , g]
      }
    }

    for (g in 1:k) {
      tau_value <- sqrt(sum(temp[, g]) / (sum(z[, g]) * d))
      tau[, , g] <- diag(tau_value, d)
      P_fin[, , g] <- P_est[, , g]
    }
  } else if (model_type == "E.comp") {
    ### P ###

    M_g <- array(NA, dim = c(d, n, k))

    for (g in 1:k) {
      Z_g[, , g] <- diag(z[, g])

      for (i in 1:n) {
        C_g[, i, g] <- diag(sqrt(w[i, , g])) %*% (X[i, ] - mu[, g])
      }

      M_g[, , g] <- solve(diag(sigma[, g])) %*% C_g[, , g]

      P_init[, g] <- P2l(P[, , g])
    }

    P_init <- as.vector(P_init)

    if (method == "L-BFGS-B") {
      opt_res <- optim(
        par = P_init, fn = obj_fun_Pg_E.comp, method = method, control = list(fnscale = -1, maxit = maxit, trace = trace),
        M = M_g, Z = Z_g, d = d, k = k, z = z, lower = lower, upper = upper
      )
    } else {
      opt_res <- optim(
        par = P_init, fn = obj_fun_Pg_E.comp, method = method, control = list(fnscale = -1, maxit = maxit, trace = trace),
        M = M_g, Z = Z_g, d = d, k = k, z = z
      )
    }

    P_opt2 <- matrix(opt_res$par, d * (d - 1) / 2, k)
    P_est <- array(NA, dim = c(d, d, k))

    for (g in 1:k) {
      P_est[, , g] <- l2P(P_opt2[, g])
      P_est[, , g] <- (P_est[, , g] + t(P_est[, , g])) / 2
    }

    ### T ###

    S <- array(NA, dim = c(d, d, n, k))

    for (g in 1:k) {
      for (i in 1:n) {
        S[, , i, g] <- z[i, g] * diag(as.vector(diag(sqrt(w[i, , g])) %*% (X[i, ] - mu[, g])), d) %*% solve(P_est[, , g]) %*% diag(as.vector(diag(sqrt(w[i, , g])) %*% (X[i, ] - mu[, g])), d)
      }
    }

    T_init <- log(sigma[, 1])

    opt_res2 <- optim(
      par = T_init, fn = obj_fun_Pg_E.comp_tau, method = method2, control = list(fnscale = -1, maxit = maxit, trace = trace),
      n = n, S = rowSums(S, dims = 2), d = d
    )

    T_opt <- exp(opt_res2$par)

    # Get the optimized parameters and transform them back into the matrix

    for (g in 1:k) {
      tau[, , g] <- diag(T_opt, d)
      P_fin[, , g] <- P_est[, , g]
    }
  } else if (model_type == "Free") {
    c <- array(NA, dim = c(n, d, k))

    for (g in 1:k) {
      for (i in 1:n) {
        c[i, , g] <- diag(sqrt(w[i, , g])) %*% (X[i, ] - mu[, g])
      }
    }

    temp <- array(NA, dim = c(d, d, n, k))

    for (i in 1:n) {
      for (g in 1:k) {
        temp[, , i, g] <- z[i, g] * c[i, , g] %*% t(c[i, , g])
      }
    }

    Sig <- array(NA, dim = c(d, d, k))
    Sig.dec <- list()

    for (g in 1:k) {
      Sig[, , g] <- rowSums(temp[, , , g], dims = 2) / sum(z[, g])
      Sig.dec[[g]] <- scale_correlation_decomp(Sig[, , g])

      tau[, , g] <- Sig.dec[[g]]$T
      P_fin[, , g] <- Sig.dec[[g]]$P
    }
  }

  return(list(P = P_fin, tau = tau))
}

#### Function to update theta ####

update_theta <- function(z, w, n, d, k, model_type, Asympt = 20) {
  theta <- matrix(NA, d, k)
  temp <- array(NA, dim = c(n, d, k))

  for (i in 1:n) {
    for (j in 1:d) {
      for (g in 1:k) {
        temp[i, j, g] <- z[i, g] * (w[i, j, g] - 1)
      }
    }
  }

  if (model_type == "E.comp.dim") {
    theta_value <- (n * d) / sum(temp)

    for (g in 1:k) {
      theta[, g] <- rep(theta_value, d)
    }
  } else if (model_type == "E.dim") {
    for (g in 1:k) {
      theta[, g] <- rep((sum(z[, g]) * d) / sum(temp[, , g]), d)
    }
  } else if (model_type == "E.comp") {
    for (j in 1:d) {
      theta[j, ] <- n / sum(temp[, j, ])
    }
  } else if (model_type == "Free") {
    for (j in 1:d) {
      for (g in 1:k) {
        theta[j, g] <- sum(z[, g]) / sum(temp[, j, g])
      }
    }
  } else if (model_type == "Asympt") {
    theta <- matrix(Asympt, d, k)
  }

  return(theta)
}

#### Density DSSEN ####

dDSSEN <- function(x, mu = rep(0, d), Sigma, theta = Inf) {
  if (missing(Sigma)) {
    stop("Sigma is missing")
  }
  if (sum(theta < 0) > 0) {
    stop("theta must be greater than, or equal to, 0")
  }

  if (is.matrix(Sigma)) {
    d <- ncol(Sigma)
    q <- nrow(x)
  }

  if (!is.matrix(Sigma)) {
    d <- 1
    q <- length(x)
  }

  if (is.vector(x)) {
    x <- matrix(x, length(x), 1)
    Sigma <- matrix(Sigma, nrow = d, ncol = d)
  }

  PDF <- NULL

  Sigmainv <- solve(Sigma)
  const1 <- 2^d / (sqrt(det(Sigma))) * prod(theta) * prod(exp(theta))

  if (d == 1) {
    Dt <- theta

    fpass <- function(z) {
      Dy <- z - mu
      P <- Dy %*% Sigmainv %*% Dy
      Pt <- P + 2 * Dt
      Ptinv <- solve(Pt)
      const2 <- sqrt(det(Ptinv))

      M <- as.vector(mnormt::recintab(rep(2, d), a = rep(1, d), b = rep(Inf, d), mu = rep(0, d), S = Ptinv))

      PDFpass <- const1 * const2 * abs(M[length(M)])

      return(PDFpass)
    }
  }

  if (d > 1) {
    Dt <- diag(theta)

    fpass <- function(z) {
      Dy <- diag(z - mu)
      P <- Dy %*% Sigmainv %*% Dy
      Pt <- P + 2 * Dt
      Ptinv <- solve(Pt)
      const2 <- sqrt(det(Ptinv))

      M <- as.vector(mnormt::recintab(rep(2, d), a = rep(1, d), b = rep(Inf, d), mu = rep(0, d), S = Ptinv))

      PDFpass <- const1 * const2 * abs(M[length(M)])

      return(PDFpass)
    }
  }

  PDF <- sapply(1:q, function(i) fpass(z = x[i, ]))

  return(PDF)
}

#### Random generation DSSEN ####

rDSSEN <- function(n, mu = rep(0, d), Sigma = NULL, theta = rep(500, d), norm.dens = "dmnorm") {
  # Lambda is a vector with diagonal elements

  if (min(theta) < 0) {
    stop("each element in theta must be greater than, or equal to, 0")
  }

  if (is.vector(mu)) {
    d <- length(mu)
  }
  if (!is.vector(mu)) {
    d <- 1
  }

  X <- w <- array(0, c(n, d), dimnames = list(1:n, paste("X.", 1:d, sep = "")))
  for (h in 1:d) {
    w[, h] <- (1 + stats::rexp(n = n, theta[h]))
  }

  DeltaRS <- SigmaMod <- array(0, c(n, d, d), dimnames = list(1:n, paste("X.", 1:d, sep = ""), paste("X.", 1:d, sep = "")))
  for (i in 1:n) {
    if (d > 1) {
      DeltaRS[i, , ] <- diag(1 / sqrt(w[i, ]))
    }
    if (d == 1) {
      DeltaRS[i, , ] <- 1 / sqrt(w[i, ])
    }
    SigmaMod[i, , ] <- DeltaRS[i, , ] %*% Sigma %*% DeltaRS[i, , ]
    if (norm.dens == "dmnorm") {
      X[i, ] <- mnormt::rmnorm(n = 1, mean = mu, varcov = SigmaMod[i, , ])
    }
    if (norm.dens == "dmvnorm") {
      X[i, ] <- mvtnorm::rmvnorm(n = 1, mean = mu, sigma = SigmaMod[i, , ])
    }
  }

  return(list(X = X, w = w))
}

#### From real vector of i-1 dimension to coordinates on the hypersphere of dimension i #####

real2sphere <- function(
    real # vector of i-1 dimension (real numbers) giving coordinates on hypersphere in dimension i
    ) {
  p <- length(real) + 1

  polar <- atan(real) + pi / 2
  polar.f <- c(polar, 0)

  y <- NULL

  if (p == 2) {
    y[1] <- cos(polar)
    y[2] <- sin(polar)
  }

  if (p > 2) {
    M1 <- matrix(rep(sin(polar.f), each = p), p, p)
    M2 <- M1 * lower.tri(M1)
    M3 <- M2 + diag(cos(polar.f)) + upper.tri(M1)
    M4 <- M3[, -p]
    y <- apply(M4, 1, prod)
  }

  return(y)
}

#### Check if the vector l of the matrix L has right dimension ####

is.wholenumber <- function(l, tol = .Machine$double.eps^0.5) abs(l - round(l)) < tol

#### From the vector l (free elements of L) of dimension d(d-1)/2 to the correlation matrix P ####

l2P <- function(l) {
  n <- length(l)
  d <- (1 + sqrt(1 + 8 * n)) / 2

  if (is.wholenumber(d) == 0) {
    print("Vector l has wrong dimension!!")
  }

  L <- c(1, rep(0, d - 1))

  seq.vec <- 1:n
  start.point <- c(0, cumsum(seq.vec))
  sub.vec <- list()
  for (j in 1:(d - 1)) {
    sub.vec[[j]] <- l[(start.point[j] + 1):(start.point[j + 1])]
  }

  for (i in 2:d) {
    L <- rbind(L, c(real2sphere(sub.vec[[i - 1]]), rep(0, d - i)))
  }

  P <- L %*% t(L)

  return(P)
}

#### Transforming a row of L which contain coordinates on hyphersphere to polar coordinates ####

rect2polar_vec <- function(v) {
  q_NAs <- sum(is.na(v))
  v_to_conv <- v[1:(length(v) - q_NAs)]
  u <- v_to_conv
  if (length(u) > 1) {
    u <- SphericalCubature::rect2polar(u)$phi
  }
  return(u)
}

#### From the correlation matrix P to the vector of free elements of L with dimension d(d-1)/2 ####

P2l <- function(P) {
  L <- t(chol(P))

  L[upper.tri(L)] <- NA

  l.polar <- unlist(apply(L, MARGIN = 1, FUN = rect2polar_vec, simplify = TRUE))[-1]

  l <- tan(l.polar - pi / 2)

  return(as.vector(l))
}
