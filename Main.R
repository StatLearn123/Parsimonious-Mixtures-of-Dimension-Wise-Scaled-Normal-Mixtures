library(foreach)
# library(snow) # must be installed
# library(doSNOW) # must be installed
# library(progress) # must be installed


#### Fit, via EM-based algorithms, of the DSNM mixture (parallel) ####

DSNM_M.fit <- function(X, # data matrix
                       k, # number of mixture components
                       corr = "all",
                       scale = "all",
                       tailedness = "all",
                       rel.tol = 0.001, # one of the 2 stopping rules (see max.iter for the other)
                       iter.max = 1000, # maximum number of iterations
                       verbose = TRUE,
                       nThreads = 1) {
  comb <- function(x, ...) {
    lapply(
      seq_along(x),
      function(i) c(x[[i]], lapply(list(...), function(y) y[[i]]))
    )
  }

  if (any(corr == "all")) {
    mod.corr <- c("Identity", "E.comp", "Free")
  } else {
    mod.corr <- corr
  }
  if (any(scale == "all")) {
    mod.scale <- c("E.comp.dim", "E.dim", "E.comp", "Free")
  } else {
    mod.scale <- scale
  }
  if (any(tailedness == "all")) {
    mod.tailedness <- c("E.comp.dim", "E.dim", "E.comp", "Free", "Asympt")
  } else {
    mod.tailedness <- tailedness
  }

  req.model <- expand.grid(mod.corr, mod.scale, mod.tailedness)
  names(req.model) <- c("corr", "scale", "tailedness")

  oper <- vector(mode = "list", length = length(k))
  time <- numeric(length(k))
  method <- "L-BFGS-B"
  method2 <- "Nelder-Mead"
  Asympt <- 100
  
  for (g in 1:length(k)) {
    if (verbose == TRUE) {
      print(paste("Fitting Parsimonious Mixtures of DSNMs", "with k =", k[g]))
    }

    ptm <- proc.time()
    
    cluster <- snow::makeCluster(nThreads, type = "SOCK")
    doSNOW::registerDoSNOW(cluster)

    varlist <- c(
      "scale_correlation_decomp", "npar.pars", "DSMN.clust.new", "update_tau_I",
      "obj_fun_P_E.comp.dim", "obj_fun_P_E.dim", "obj_fun_P_Free", "obj_fun_P_Free_tau",
      "update_tau_P", "obj_fun_Pg_E.comp.dim", "obj_fun_Pg_E.dim", "obj_fun_Pg_E.comp",
      "obj_fun_Pg_E.comp_tau", "update_tau_Pg", "update_theta", "dDSSEN", "real2sphere",
      "is.wholenumber", "l2P","rect2polar_vec", "P2l"
    )

    pb <- progress::progress_bar$new(
      format = "(:spin) [:bar] :percent [Elapsed time: :elapsedfull || Estimated time remaining: :eta]",
      total = nrow(req.model),
      complete = "=", # Completion bar character
      incomplete = "-", # Incomplete bar character
      current = ">", # Current bar character
      width = 100
    )
    progress <- function(n) {
      pb$tick()
    }
    opts <- list(progress = progress)

    l <- 0

    oper[[g]] <- foreach::foreach(
      l = 1:nrow(req.model), .combine = "comb", .multicombine = TRUE, .export = varlist,
      .init = list(list()), .packages = "mclust", .options.snow = opts
    ) %dopar% {
      res <- tryCatch(DSMN.clust.new(
        X = X, k = k[g], corr = req.model[l, 1], scale = req.model[l, 2], tailedness = req.model[l, 3], rel.tol = rel.tol, iter.max = iter.max, method = method, method2 = method2, Asympt = Asympt
      ), error = function(e) {
        NA
      })

      list(res)
    }

    snow::stopCluster(cluster)
    foreach::registerDoSEQ()
    
    ptm2 <- proc.time() - ptm
    time[g] <- ptm2[3]
  }

  return(list(
    results = oper,
    k = k,
    req.model = req.model,
    comp.time = time
  ))
} 

#### Select best fitting model ####

extract.bestM <- function(results, sp.th = 0, criterion = "BIC") {
  top <- 1
  k <- length(results[["results"]])
  num.mod <- length(results[["results"]][[1]][[1]])
  list.mod <- results[["req.model"]]
  list.mod2 <- do.call("rbind", replicate(k, list.mod, simplify = FALSE))
  count.k <- sort(rep(1:k, num.mod))
  count.mod <- rep(1:num.mod, k)
  list.mod3 <- data.frame(list.mod2, count.k, count.mod)
  
  all_INF <- numeric(k * num.mod)
  
  cont <- 0
  
  for (j in 1:k) {
    for (i in 1:num.mod) {
      if (!all(is.na(results[["results"]][[j]][[1]][[i]]))) {
        cont <- cont + 1
        all_INF[cont] <- -results[["results"]][[j]][[1]][[i]][[criterion]]
        
        kk <- results[["results"]][[j]][[1]][[i]][["prior"]]
        
        if (any(kk < sp.th)) {
          all_INF[cont] <- NA
        }
      } else {
        cont <- cont + 1
        all_INF[cont] <- NA
      }
    }
  }
  
  topBIC <- which(all_INF >= sort(all_INF, decreasing = T)[top])
  topBIC.order <- order(all_INF[topBIC], decreasing = T)
  tempBIC <- list.mod3[topBIC[topBIC.order], ]
  bestBIC <- vector(mode = "list", length = top)
  for (i in 1:top) {
    bestBIC[[i]] <- results[["results"]][[as.numeric(tempBIC[i, 4])]][[1]][[as.numeric(tempBIC[i, 5])]]
  }
  
  return(best = bestBIC)
} # extract best fitting model