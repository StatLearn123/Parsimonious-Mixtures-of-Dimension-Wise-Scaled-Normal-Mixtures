### Simulation 2 ###

library(ContaminatedMixt)
library(teigen)
library(leptokurticMixture)
library(mixSPE)
library(SenTinMixt)
library(pgmm)
library(tcltk)

n  <- 100          # obs per group
d  <- 7            # separation of the centers
df <- 3            # degrees of freedom

mu1 <- c(-d/2, 0)
mu2 <- c( d/2, 0)

nsims <- 100

## Generating Function ##

gen_cluster_heavytail <- function(n, mu, df) {
  
  U <- rt(n, df = df)      
  V <- rt(n, df = df)
  
  X1 <- U
  X2 <- tanh(U) + 0.6 * V
  
  X <- cbind(X1, X2)
  
  X <- sweep(X, 2, mu, "+")
  
  return(X)
}

pb <- tkProgressBar(
  title = "In Progress",
  label = "Start...",
  min = 0,
  max = nsims,
  initial = 0
)

#### Fitting  ####

X <- group <- results <- resMN <- resMT <- resCN <- resLN <- resPE <- resSEN <- resTIN <- resMN.F <- vector(mode = "list", length = nsims)
nThreads <- 28 # for parallel computing of DSEN_M.fit

for (l in 1:nsims) {
  
  setTkProgressBar(pb, l, label = paste("Fitting dataset number", l))
  
  ##  Data Generation
  
  X1 <- gen_cluster_heavytail(n, mu1, df)
  X2 <- gen_cluster_heavytail(n, mu2, df)
  
  X[[l]] <- rbind(X1, X2)
  group[[l]] <- factor(rep(c("G1", "G2"), each = n))
  
  #### Multivariate DSSEN Mixtures (Eigen) ####
  
  results[[l]] <- tryCatch(DSNM_M.fit(X=X[[l]], k=2, corr = "all", scale = "all", tailedness = "all", nThreads = nThreads, verbose = T), error = function(e) {NA})
  
  #### Multivariate Normal Mixtures (Eigen) ####
  
  resMN[[l]] <- vector(mode = "list", length = 4) 
  
  set.seed(123)
  resMN[[l]][[1]]<- tryCatch(CNmixt(X=X[[l]], G=2, contamination = FALSE, initialization = "mixt"), error = function(e) {NA})
  set.seed(123)
  resMN[[l]][[2]]<- tryCatch(CNmixt(X=X[[l]], G=2, contamination = FALSE, initialization = "kmeans"), error = function(e) {NA})
  set.seed(123)
  resMN[[l]][[3]]<- tryCatch(CNmixt(X=X[[l]], G=2, contamination = FALSE, initialization = "random.post"), error = function(e) {NA})
  set.seed(123)
  resMN[[l]][[4]]<- tryCatch(CNmixt(X=X[[l]], G=2, contamination = FALSE, initialization = "random.clas"), error = function(e) {NA})
  
  #### Multivariate t Mixtures (Eigen) ####
  
  resMT[[l]] <- vector(mode = "list", length = 4) 
  
  set.seed(123)
  resMT[[l]][[1]]<- tryCatch(teigen(x=X[[l]], Gs=2, scale = FALSE, dfupdate = "numeric", init = "kmeans"), error = function(e) {NA})
  set.seed(123)
  resMT[[l]][[2]]<- tryCatch(teigen(x=X[[l]], Gs=2, scale = FALSE, dfupdate = "numeric", init = "hard"), error = function(e) {NA})
  set.seed(123)
  resMT[[l]][[3]]<- tryCatch(teigen(x=X[[l]], Gs=2, scale = FALSE, dfupdate = "numeric", init = "soft"), error = function(e) {NA})
  set.seed(123)
  resMT[[l]][[4]]<- tryCatch(teigen(x=X[[l]], Gs=2, scale = FALSE, dfupdate = "numeric", init = "emem", ememargs = list(numstart=25,iter=5,model="UUUU",init="hard")), error = function(e) {NA})
  
  #### Multivariate CN Mixtures (Eigen) ####
  
  resCN[[l]] <- vector(mode = "list", length = 4) 
  
  set.seed(123)
  resCN[[l]][[1]]<- tryCatch(CNmixt(X=X[[l]], G=2, contamination = TRUE, initialization = "mixt"), error = function(e) {NA})
  set.seed(123)
  resCN[[l]][[2]]<- tryCatch(CNmixt(X=X[[l]], G=2, contamination = TRUE, initialization = "kmeans"), error = function(e) {NA})
  set.seed(123)
  resCN[[l]][[3]]<- tryCatch(CNmixt(X=X[[l]], G=2, contamination = TRUE, initialization = "random.post"), error = function(e) {NA})
  set.seed(123)
  resCN[[l]][[4]]<- tryCatch(CNmixt(X=X[[l]], G=2, contamination = TRUE, initialization = "random.clas"), error = function(e) {NA})
  
  #### Multivariate LN Mixtures (Eigen) ####
  
  resLN[[l]] <-vector(mode = "list", length = 28) 
  
  covmodels <- c("EII", "VII", "EEI", "VEI", "EVI" ,"VVI", "EEE" ,"EEV" ,"VEV", "VVV" ,"EVE", "VVE", "VEE" ,"EVV") 
  
  for (i in 1:14) {
    
    set.seed(123)
    resLN[[l]][[i]] <- tryCatch(pmln(data = X[[l]], G = 2, kml = c(1,1,1), scale.data = FALSE, pprogress = TRUE, covModels = covmodels[i], betaModels = "V"), error = function(e) {NA})
    
  }
  for (i in 15:28) {
    
    set.seed(123)
    resLN[[l]][[i]] <- tryCatch(pmln(data = X[[l]], G = 2, kml = c(1,1,1), scale.data = FALSE, pprogress = TRUE, covModels = covmodels[i-14], betaModels = "E"), error = function(e) {NA})
    
  }
  
  #### Multivariate PE Mixtures (Eigen) ####
  
  set.seed(123)
  resPE[[l]] <- tryCatch(EMGr(data = X[[l]], G = 2, scale = FALSE, keepResults = TRUE, initialization = 1), error = function(e) {NA})
  
  #### Multivariate SEN Mixtures (Eigen) ####
  
  set.seed(123)
  resSEN.init <- tryCatch(Mixt.fit.init(X=X[[l]], k=2, density = "MSEN", ncores = 10, verbose = TRUE, nstartR = 50), error = function(e) {NA})
  resSEN[[l]] <- tryCatch(Mixt.fit(X=X[[l]], k = 2, init.par = resSEN.init, density = "MSEN", ncores = 28, verbose = T), error = function(e) {NA})
  
  #### Multivariate TIN Mixtures (Eigen) ####
  
  set.seed(123)
  resTIN.init <- tryCatch(Mixt.fit.init(X=X[[l]], k=2, density = "MTIN", ncores = 10, verbose = TRUE, nstartR = 50), error = function(e) {NA})
  resTIN[[l]] <- tryCatch(Mixt.fit(X=X[[l]], k = 2, init.par = resTIN.init, density = "MTIN", ncores = 28, verbose = T), error = function(e) {NA})
  
  #### Multivariate Normal Mixtures (Factor) ####
  
  resMN.F[[l]] <- vector(mode = "list", length = 4) 
  resMN.F[[l]][[1]] <- tryCatch(pgmmEM(x=X[[l]], rG=2, rq=1, zstart = 1, relax = TRUE), error = function(e) {NA})
  resMN.F[[l]][[2]] <- tryCatch(pgmmEM(x=X[[l]], rG=2, rq=1, zstart = 2, relax = TRUE), error = function(e) {NA})
  
  resMN.F[[l]][[3]] <- tryCatch(pgmmEM(x=X[[l]], rG=2, rq=1, zstart = 1, relax = TRUE, icl = TRUE), error = function(e) {NA})
  resMN.F[[l]][[4]] <- tryCatch(pgmmEM(x=X[[l]], rG=2, rq=1, zstart = 2, relax = TRUE, icl = TRUE), error = function(e) {NA})
  
  rm(.Random.seed, envir=globalenv())

}

close(pb)

#### Extract the results ####

library(mclust)

res_fit <- res_ARI <- array(NA, dim = c(9,2,nsims), dimnames = list(c("DSSEN","MN","MT","MCN","MLN","MPE","MSEN","MTIN","MN-F"),
                                                                    c("BIC","ICL"),
                                                                    1:nsims)) # 9 models, 2 ic, nsims datasets
k <- length(unique(group[[1]]))
num <- n*k

for (l in 1:nsims) {
  
  trueclass<- group[[l]]
  
  ### Multivariate DSSEN Mixtures (Eigen) ###
  
  win     <- extract.bestM(results = results[[l]], criterion = "BIC")
  win2    <- extract.bestM(results = results[[l]], criterion = "ICL")
  
  res_fit[1,1,l] <- win[[1]][["BIC"]]
  res_fit[1,2,l] <- win[[1]][["ICL"]]
  res_ARI[1,1,l] <- round(adjustedRandIndex(trueclass,win[[1]][["clust.vec"]]), digits = 2)
  res_ARI[1,2,l] <- round(adjustedRandIndex(trueclass,win[[1]][["clust.vec"]]), digits = 2)
  
  ### Multivariate Normal Mixtures (Eigen) ###
  
  resMNbic <- resMNicl <- numeric(4)
  winMN.BIC <- winMN.ICL <-vector(mode = "list", length = 4) 
  
  for (i in 1:4) {
    
    winMN.BIC[[i]] <- tryCatch(getBestModel(resMN[[l]][[i]], criterion = "BIC"), error = function(e) {NA})
    winMN.ICL[[i]] <- tryCatch(getBestModel(resMN[[l]][[i]], criterion = "ICL"), error = function(e) {NA})
    
    resMNbic[i] <- tryCatch(winMN.BIC[[i]][["models"]][[1]][["IC"]][["BIC"]], error = function(e) {NA})
    resMNicl[i] <- tryCatch(winMN.ICL[[i]][["models"]][[1]][["IC"]][["ICL"]], error = function(e) {NA})
    
  }
  
  res_fit[2,1,l] <- -resMNbic[[which.max(resMNbic)]]
  res_fit[2,2,l] <- -resMNicl[[which.max(resMNicl)]]
  res_ARI[2,1,l] <- round(adjustedRandIndex(trueclass,winMN.BIC[[which.max(resMNbic)]][["models"]][[1]][["group"]]), digits = 2)
  res_ARI[2,2,l] <- round(adjustedRandIndex(trueclass,winMN.ICL[[which.max(resMNicl)]][["models"]][[1]][["group"]]), digits = 2)
  
  #### Multivariate t Mixtures (Eigen) ####
  
  resMTbic <- resMTicl <- numeric(4)
  
  for (i in 1:4) {
    
    resMTbic[i] <- resMT[[l]][[i]][["bic"]]
    resMTicl[i] <- resMT[[l]][[i]][["iclresults"]][["icl"]]
    
  }
  
  res_fit[3,1,l] <- -resMT[[l]][[which.max(resMTbic)]][["bic"]] 
  res_fit[3,2,l] <- -resMT[[l]][[which.max(resMTicl)]][["iclresults"]][["icl"]]
  res_ARI[3,1,l] <- round(adjustedRandIndex(trueclass,resMT[[l]][[which.max(resMTbic)]][["classification"]] ), digits = 2)
  res_ARI[3,2,l] <- round(adjustedRandIndex(trueclass,resMT[[l]][[which.max(resMTicl)]][["iclresults"]][["classification"]]), digits = 2)
  
  ### Multivariate CN Mixtures (Eigen) ###
  
  resCNbic <- resCNicl <- numeric(4)
  winCN.BIC <- winCN.ICL <-vector(mode = "list", length = 4) 
  
  for (i in 1:4) {
    
    winCN.BIC[[i]] <- tryCatch(getBestModel(resCN[[l]][[i]], criterion = "BIC"), error = function(e) {NA})
    winCN.ICL[[i]] <- tryCatch(getBestModel(resCN[[l]][[i]], criterion = "ICL"), error = function(e) {NA})
    
    resCNbic[i] <- tryCatch(winCN.BIC[[i]][["models"]][[1]][["IC"]][["BIC"]], error = function(e) {NA})
    resCNicl[i] <- tryCatch(winCN.ICL[[i]][["models"]][[1]][["IC"]][["ICL"]], error = function(e) {NA})
    
  }
  
  res_fit[4,1,l] <- -resCNbic[[which.max(resCNbic)]]
  res_fit[4,2,l] <- -resCNicl[[which.max(resCNicl)]]
  res_ARI[4,1,l] <- round(adjustedRandIndex(trueclass,winCN.BIC[[which.max(resCNbic)]][["models"]][[1]][["group"]]), digits = 2)
  res_ARI[4,2,l] <- round(adjustedRandIndex(trueclass,winCN.ICL[[which.max(resCNicl)]][["models"]][[1]][["group"]]), digits = 2)
  
  ### Multivariate LN Mixtures (Eigen) ###
  
  resLNbic <- resLNicl <- numeric(28)
  for (i in 1:28) {
    
    resLNbic[i] <- tryCatch(resLN[[l]][[i]][["bicModel"]][["bic"]], error = function(e) {NA})
    
    z <- tryCatch(t(resLN[[l]][[i]][["z"]]), error = function(e) {NA})
    z.const <- tryCatch((z < 10^(-322)) * 10^(-322) + (z > 10^(-322)) * z , error = function(e) {NA}) 
    hard.z <- tryCatch((matrix(rep(apply(z, 1, max), k), num, k, byrow = F) == z) * 1, error = function(e) {NA})
    EN <- tryCatch(-sum(hard.z * log(z.const)), error = function(e) {NA})
    
    resLNicl[i] <- tryCatch(resLNbic[i] + 2 * EN, error = function(e) {NA})
    
  }
  
  res_fit[5,1,l] <- resLNbic[[which.min(resLNbic)]]
  res_fit[5,2,l] <- resLNicl[[which.min(resLNicl)]]
  res_ARI[5,1,l] <- round(adjustedRandIndex(trueclass,resLN[[l]][[which.min(resLNbic)]][["map"]]), digits = 2)
  res_ARI[5,2,l] <- round(adjustedRandIndex(trueclass,resLN[[l]][[which.min(resLNicl)]][["map"]]), digits = 2)
  
  ### Multivariate PE Mixtures (Eigen) ###
  
  resPEbic <- resPEicl <- numeric(16)
  for (i in 1:16) {
    
    resPEbic[i] <- resPE[[l]][["BIC"]][i]
    
    z <- tryCatch(resPE[[l]][["allModels"]][[k]][[i]][["z"]], error = function(e) {NA})
    
    z.const <- tryCatch((z < 10^(-322)) * 10^(-322) + (z > 10^(-322)) * z , error = function(e) {NA}) 
    hard.z <- tryCatch((matrix(rep(apply(z, 1, max), k), num, k, byrow = F) == z) * 1, error = function(e) {NA})
    EN <- tryCatch(-sum(hard.z * log(z.const)), error = function(e) {NA})
    
    resPEicl[i] <- tryCatch(resPEbic[i] + 2 * EN, error = function(e) {NA})
    
  }
  
  res_fit[6,1,l] <- resPEbic[[which.min(resPEbic)]]
  res_fit[6,2,l] <- resPEicl[[which.min(resPEicl)]]
  res_ARI[6,1,l] <- round(adjustedRandIndex(trueclass,resPE[[l]][["allModels"]][[k]][[which.min(resPEbic)]][["map"]]), digits = 2)
  res_ARI[6,2,l] <- round(adjustedRandIndex(trueclass,resPE[[l]][["allModels"]][[k]][[which.min(resPEbic)]][["map"]]), digits = 2)
  
  ### Multivariate SEN Mixtures (Eigen) ###
  
  res_fit[7,1,l] <- resSEN[[l]][["BicWin"]][["BIC"]]
  res_fit[7,2,l] <- resSEN[[l]][["IclWin"]][["ICL"]]
  res_ARI[7,1,l] <- round(adjustedRandIndex(trueclass,resSEN[[l]][["BicWin"]][["classification"]]), digits = 2)
  res_ARI[7,2,l] <- round(adjustedRandIndex(trueclass,resSEN[[l]][["IclWin"]][["classification"]]), digits = 2)
  
  ### Multivariate TIN Mixtures (Eigen) ###
  
  res_fit[8,1,l] <- resTIN[[l]][["BicWin"]][["BIC"]]
  res_fit[8,2,l] <- resTIN[[l]][["IclWin"]][["ICL"]]
  res_ARI[8,1,l] <- round(adjustedRandIndex(trueclass,resTIN[[l]][["BicWin"]][["classification"]]), digits = 2)
  res_ARI[8,2,l] <- round(adjustedRandIndex(trueclass,resTIN[[l]][["IclWin"]][["classification"]]), digits = 2)
  
  ### Multivariate Normal Mixtures (Factor) ###
  
  res_fit[9,1,l] <- -max(c(resMN.F[[l]][[1]][["summ_info"]][[4]],resMN.F[[l]][[2]][["summ_info"]][[4]]))
  res_fit[9,2,l] <- -max(c(resMN.F[[l]][[3]][["summ_info"]][[4]],resMN.F[[l]][[4]][["summ_info"]][[4]]))
  res_ARI[9,1,l] <- round(adjustedRandIndex(trueclass,resMN.F[[l]][[which.max(max(c(resMN.F[[l]][[1]][["summ_info"]][[4]],resMN.F[[l]][[2]][["summ_info"]][[4]])))]][["map"]]), digits = 2)
  res_ARI[9,2,l] <- round(adjustedRandIndex(trueclass,resMN.F[[l]][[which.max(max(c(resMN.F[[l]][[3]][["summ_info"]][[4]],resMN.F[[l]][[4]][["summ_info"]][[4]])))]][["map"]]), digits = 2)
  
}

## BIC, ICL e ARI averages ##

apply(res_fit, c(1, 2), mean)
apply(res_ARI, c(1, 2), mean)




