### Simulation 1 - Example ###

n <- 200  # Number of observations
d <- 3    # Number of dimensions
k <- 2    # Number of clusters

mu <- matrix(NA, d, k) # Cluster means (d x k matrix)
Sigma <- tau <- array(NA, dim = c(d,d,k))
theta <- matrix(NA, d, k)

## I / E.C.D / E.C.D. ##

P <- array(diag(d), dim = c(d,d,k))
tau[,,1] <- tau[,,2] <- diag(rep(1.10,d), d) 
theta <- matrix(0.1, d, k)

mu[,1] <- c(1,2,1)
mu[,2] <- mu[,1] + 2

Sigma[,,1] <- tau[,,1] %*% P[,,1] %*% tau[,,1]
Sigma[,,2] <- tau[,,2] %*% P[,,2] %*% tau[,,2]

pii <- rep(1/k,k)

### Data Generation ###

n_iter <- 1
nThreads = 28 # number of cores per parallel computing
res <- X <- sizes <- vector(mode = "list", length = n_iter)

## Load the parameters here ##

for (j in 1:n_iter) {
  
  print(paste("Iteration number", j))
  w <- rmultinom(n, size = 1, prob = pii) 
  sizes[[j]] <- rowSums(w)
  
  X1 <- rDSSEN(n=sizes[[j]][1], mu=mu[,1], Sigma = Sigma[,,1], theta = theta[,1]) 
  X2 <- rDSSEN(n=sizes[[j]][2], mu=mu[,2], Sigma = Sigma[,,2], theta = theta[,2]) 
  
  X[[j]] <- rbind(X1$X, X2$X)
  
  res[[j]] <- tryCatch(DSNM_M.fit(X=X[[j]], k=k, corr = "all", scale = "all", tailedness = "all", nThreads = nThreads, verbose = F), error = function(e) {NA}) 
  # win <- tryCatch(extract.bestM(results = res[[j]]), error = function(e) {NA})
  
}


