library(Stat2Data)
library(dplyr)
library(tidyr)

data(Hawks)

data1 <- select(Hawks, Species, Sex, Wing, Weight, Culmen, Tail)
data2 <- filter(data1, Species %in% c("SS","CH"))
data3 <- filter(data2, Sex %in% c("F","M"))
data4 <- data3[,c(3:6,1:2)]
variab<- c(1:4)

X  <- data4[,variab]
X  <- X %>% drop_na()
X  <- as.matrix(X)

nThreads <- 28 # for parallel computing of DSEN_M.fit

results <- DSNM_M.fit(X=X, k=4, corr = "all", scale = "all", tailedness = "all", nThreads = nThreads, verbose = T)
win     <- extract.bestM(results = results, criterion = "BIC")
win2     <- extract.bestM(results = results, criterion = "ICL")

### Classification ###

library(mclust)

Y  <- data4
Y  <- Y %>% drop_na()
trueclass <- numeric(length = nrow(Y))
for (i in 1:nrow(Y)) {
  
  if(Y$Sex[i]=="M" & Y$Species[i]=="CH"){
    trueclass[i] <- 1
  }
  if(Y$Sex[i]=="M" & Y$Species[i]=="SS"){
    
    trueclass[i] <- 2
  }
  if(Y$Sex[i]=="F" & Y$Species[i]=="CH"){
    trueclass[i] <- 3
  }
  if(Y$Sex[i]=="F" & Y$Species[i]=="SS"){
    
    trueclass[i] <- 4
  }
  
}

round(mclust::adjustedRandIndex(trueclass, win[[1]]$clust.vec), digits = 2)

### Alternative models ###

library(ContaminatedMixt)
library(teigen)
library(leptokurticMixture)
library(mixSPE)
library(SenTinMixt)
library(pgmm)

resAlt <- resAlt.ARI <- matrix(NA, nrow = 8, ncol = 2)
rownames(resAlt) <- rownames(resAlt.ARI) <- c("MN","MT","MCN","MLN","MPE","MSEN","MTIN","MN-F")
colnames(resAlt) <- colnames(resAlt.ARI) <- c("BIC","ICL")

#### Multivariate Normal Mixtures (Eigen) ####

resMN <- winMN.BIC <- winMN.ICL <-vector(mode = "list", length = 4) 

set.seed(123)
resMN[[1]]<- CNmixt(X=X, G=4, contamination = FALSE, initialization = "mixt")
set.seed(123)
resMN[[2]]<- CNmixt(X=X, G=4, contamination = FALSE, initialization = "kmeans")
set.seed(123)
resMN[[3]]<- CNmixt(X=X, G=4, contamination = FALSE, initialization = "random.post")
set.seed(123)
resMN[[4]]<- CNmixt(X=X, G=4, contamination = FALSE, initialization = "random.clas")

resMNbic <- resMNicl <- numeric(4)

for (i in 1:4) {
  
  winMN.BIC[[i]] <- getBestModel(resMN[[i]], criterion = "BIC")
  winMN.ICL[[i]] <- getBestModel(resMN[[i]], criterion = "ICL")
  
  resMNbic[i] <- winMN.BIC[[i]][["models"]][[1]][["IC"]][["BIC"]]
  resMNicl[i] <- winMN.ICL[[i]][["models"]][[1]][["IC"]][["ICL"]]
  
}

resAlt[1,1] <- -resMNbic[[which.max(resMNbic)]]
resAlt[1,2] <- -resMNicl[[which.max(resMNicl)]]

resAlt.ARI[1,1] <- round(adjustedRandIndex(trueclass,winMN.BIC[[which.max(resMNbic)]][["models"]][[1]][["group"]]), digits = 2)
resAlt.ARI[1,2] <- round(adjustedRandIndex(trueclass,winMN.ICL[[which.max(resMNicl)]][["models"]][[1]][["group"]]), digits = 2)

#### Multivariate t Mixtures (Eigen) ####

resMT <- vector(mode = "list", length = 4) 

set.seed(123)
resMT[[1]]<- teigen(x=X, Gs=4, scale = FALSE, dfupdate = "numeric", init = "kmeans")
set.seed(123)
resMT[[2]]<- teigen(x=X, Gs=4, scale = FALSE, dfupdate = "numeric", init = "hard")
set.seed(123)
resMT[[3]]<- teigen(x=X, Gs=4, scale = FALSE, dfupdate = "numeric", init = "soft")
set.seed(123)
resMT[[4]]<- teigen(x=X, Gs=4, scale = FALSE, dfupdate = "numeric", init = "emem", ememargs = list(numstart=25,iter=5,model="UUUU",init="hard"))

resMTbic <- resMTicl <- numeric(4)
for (i in 1:4) {
  
  resMTbic[i] <- resMT[[i]][["bic"]]
  resMTicl[i] <- resMT[[i]][["iclresults"]][["icl"]]
  
}

resAlt[2,1] <- -resMT[[which.max(resMTbic)]][["bic"]] 
resAlt[2,2] <- -resMT[[which.max(resMTicl)]][["iclresults"]][["icl"]]

resAlt.ARI[2,1] <- round(adjustedRandIndex(trueclass,resMT[[which.max(resMTbic)]][["classification"]] ), digits = 2)
resAlt.ARI[2,2] <- round(adjustedRandIndex(trueclass,resMT[[which.max(resMTicl)]][["iclresults"]][["classification"]]), digits = 2)

#### Multivariate CN Mixtures (Eigen) ####

resCN <- winCN.BIC <- winCN.ICL <-vector(mode = "list", length = 4) 

set.seed(123)
resCN[[1]]<- CNmixt(X=X, G=4, contamination = TRUE, initialization = "mixt")
set.seed(123)
resCN[[2]]<- CNmixt(X=X, G=4, contamination = TRUE, initialization = "kmeans")
set.seed(123)
resCN[[3]]<- CNmixt(X=X, G=4, contamination = TRUE, initialization = "random.post")
set.seed(123)
resCN[[4]]<- CNmixt(X=X, G=4, contamination = TRUE, initialization = "random.clas")

resCNbic <- resCNicl <- numeric(4)

for (i in 1:4) {
  
  winCN.BIC[[i]] <- getBestModel(resCN[[i]], criterion = "BIC")
  winCN.ICL[[i]] <- getBestModel(resCN[[i]], criterion = "ICL")
  
  resCNbic[i] <- winCN.BIC[[i]][["models"]][[1]][["IC"]][["BIC"]]
  resCNicl[i] <- winCN.ICL[[i]][["models"]][[1]][["IC"]][["ICL"]]
  
}

resAlt[3,1] <- -resCNbic[[which.max(resCNbic)]]
resAlt[3,2] <- -resCNicl[[which.max(resCNicl)]]

resAlt.ARI[3,1] <- round(adjustedRandIndex(trueclass,winCN.BIC[[which.max(resCNbic)]][["models"]][[1]][["group"]]), digits = 2)
resAlt.ARI[3,2] <- round(adjustedRandIndex(trueclass,winCN.ICL[[which.max(resCNicl)]][["models"]][[1]][["group"]]), digits = 2)

#### Multivariate LN Mixtures (Eigen) ####

resLN <-vector(mode = "list", length = 28) 

covmodels <- c("EII", "VII", "EEI", "VEI", "EVI" ,"VVI", "EEE" ,"EEV" ,"VEV", "VVV" ,"EVE", "VVE", "VEE" ,"EVV") 

for (i in 1:14) {
  
  set.seed(123)
  resLN[[i]] <- pmln(data = X, G = 4, kml = c(1,1,1), scale.data = FALSE, pprogress = TRUE, covModels = covmodels[i], betaModels = "V")
  
}
for (i in 15:28) {
  
  set.seed(123)
  resLN[[i]] <- pmln(data = X, G = 4, kml = c(1,1,1), scale.data = FALSE, pprogress = TRUE, covModels = covmodels[i-14], betaModels = "E")
  
}

resLNbic <- resLNicl <- numeric(28)
for (i in 1:28) {
  
  resLNbic[i] <- tryCatch(resLN[[i]][["bicModel"]][["bic"]], error = function(e) {NA})
  
  z <- tryCatch(t(resLN[[i]][["z"]]), error = function(e) {NA})
  z.const <- tryCatch((z < 10^(-322)) * 10^(-322) + (z > 10^(-322)) * z , error = function(e) {NA})
  hard.z <- tryCatch((matrix(rep(apply(z, 1, max), 4), 322, 4, byrow = F) == z) * 1, error = function(e) {NA})
  EN <- tryCatch(-sum(hard.z * log(z.const)), error = function(e) {NA})
  
  resLNicl[i] <- tryCatch(resLNbic[i] + 2 * EN, error = function(e) {NA})
  
}

resAlt[4,1] <- resLNbic[[which.min(resLNbic)]]
resAlt[4,2] <- resLNicl[[which.min(resLNicl)]]

resAlt.ARI[4,1] <- round(adjustedRandIndex(trueclass,resLN[[which.min(resLNbic)]][["map"]]), digits = 2)
resAlt.ARI[4,2] <- round(adjustedRandIndex(trueclass,resLN[[which.min(resLNicl)]][["map"]]), digits = 2)

#### Multivariate PE Mixtures (Eigen) ####

set.seed(123)
resPE <- EMGr(data = X, G = 4, scale = FALSE, keepResults = TRUE, initialization = 1)

resPEbic <- resPEicl <- numeric(16)
for (i in 1:16) {
  
  resPEbic[i] <- resPE[["BIC"]][i]
  
  z <- tryCatch(resPE[["allModels"]][[4]][[i]][["z"]], error = function(e) {NA})
  
  z.const <- tryCatch((z < 10^(-322)) * 10^(-322) + (z > 10^(-322)) * z , error = function(e) {NA})
  hard.z <- tryCatch((matrix(rep(apply(z, 1, max), 4), 322, 4, byrow = F) == z) * 1, error = function(e) {NA})
  EN <- tryCatch(-sum(hard.z * log(z.const)), error = function(e) {NA})
  
  resPEicl[i] <- tryCatch(resPEbic[i] + 2 * EN, error = function(e) {NA})
  
}

resAlt[5,1] <- resPEbic[[which.min(resPEbic)]]
resAlt[5,2] <- resPEicl[[which.min(resPEicl)]]

resAlt.ARI[5,1] <- round(adjustedRandIndex(trueclass,resPE[["allModels"]][[4]][[which.min(resPEbic)]][["map"]]), digits = 2)
resAlt.ARI[5,2] <- round(adjustedRandIndex(trueclass,resPE[["allModels"]][[4]][[which.min(resPEbic)]][["map"]]), digits = 2)

#### Multivariate SEN Mixtures (Eigen) ####

set.seed(123)
resSEN.init <- Mixt.fit.init(X=X, k=4, density = "MSEN", ncores = 10, verbose = TRUE, nstartR = 50)
resSEN <- Mixt.fit(X=X, k = 4, init.par = resSEN.init, density = "MSEN", ncores = 28, verbose = T)

resAlt[6,1] <- resSEN[["BicWin"]][["BIC"]]
resAlt[6,2] <- resSEN[["IclWin"]][["ICL"]]

resAlt.ARI[6,1] <- round(adjustedRandIndex(trueclass,resSEN[["BicWin"]][["classification"]]), digits = 2)
resAlt.ARI[6,2] <- round(adjustedRandIndex(trueclass,resSEN[["IclWin"]][["classification"]]), digits = 2)

#### Multivariate TIN Mixtures (Eigen) ####

set.seed(123)
resTIN.init <- Mixt.fit.init(X=X, k=4, density = "MTIN", ncores = 10, verbose = TRUE, nstartR = 50)
resTIN <- Mixt.fit(X=X, k = 4, init.par = resTIN.init, density = "MTIN", ncores = 28, verbose = T)

resAlt[7,1] <- resTIN[["BicWin"]][["BIC"]]
resAlt[7,2] <- resTIN[["IclWin"]][["ICL"]]

resAlt.ARI[7,1] <- round(adjustedRandIndex(trueclass,resTIN[["BicWin"]][["classification"]]), digits = 2)
resAlt.ARI[7,2] <- round(adjustedRandIndex(trueclass,resTIN[["IclWin"]][["classification"]]), digits = 2)

#### Multivariate Normal Mixtures (Factor) ####

resMN.F <- vector(mode = "list", length = 4) 
resMN.F[[1]] <- pgmmEM(x=X, rG=4, rq=1:3, zstart = 1, relax = TRUE)
resMN.F[[2]] <- pgmmEM(x=X, rG=4, rq=1:3, zstart = 2, relax = TRUE)

resMN.F[[3]] <- pgmmEM(x=X, rG=4, rq=1:3, zstart = 1, relax = TRUE, icl = TRUE)
resMN.F[[4]] <- pgmmEM(x=X, rG=4, rq=1:3, zstart = 2, relax = TRUE, icl = TRUE)

resAlt[8,1] <- -max(c(resMN.F[[1]][["summ_info"]][[4]],resMN.F[[2]][["summ_info"]][[4]]))
resAlt[8,2] <- -max(c(resMN.F[[3]][["summ_info"]][[4]],resMN.F[[4]][["summ_info"]][[4]]))

resAlt.ARI[8,1] <- round(adjustedRandIndex(trueclass,resMN.F[[which.max(max(c(resMN.F[[1]][["summ_info"]][[4]],resMN.F[[2]][["summ_info"]][[4]])))]][["map"]]), digits = 2)
resAlt.ARI[8,2] <- round(adjustedRandIndex(trueclass,resMN.F[[which.max(max(c(resMN.F[[3]][["summ_info"]][[4]],resMN.F[[4]][["summ_info"]][[4]])))]][["map"]]), digits = 2)
