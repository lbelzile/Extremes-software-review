
## ----internal_functions, echo = FALSE------------------
mydetach <- function(j, pkglist) {
  pkglist <- pkglist[j,"package"]
  if (length(pkglist) > 0)
    lapply(
      paste0("package:", pkglist),
      # Unload add-on packages
      detach,
      character.only = TRUE,
      unload = TRUE
    )
}
myfulldetach <- function() {
  pkglist <- names(sessionInfo()$otherPkgs)
  if (length(pkglist) > 0)
    lapply(
      paste0("package:", pkglist),
      # Unload add-on packages
      detach,
      character.only = TRUE,
      unload = TRUE
    )
}
myfullunload <- function() {
  nmslist <- names(sessionInfo()$loadedOnly)
  i <- 0
  while (length(nmslist) > 0 && i < 10)
  {
    lapply(nmslist, function(x)
      try(unloadNamespace(x), silent = TRUE)
    )
    nmslist <- names(sessionInfo()$loadedOnly)
    i <- i  +  1
  }
}
# Function to use function without loading package
# https://stackoverflow.com/questions/38983179/do-call-a-function-in-r-without-loading-the-package
getfun <- function(x) {
   if(length(grep("::", x))>0) {
      parts<-strsplit(x, "::")[[1]]
      getExportedValue(parts[1], parts[2])
   } else {
      x
   }
}

check_gev_fit <- function(j,
                          obs,
                          echo = FALSE,
                          pkgfunlist,
                          pkgotherpar) {
  if(is.logical(echo)) {
    echo <- as.integer(TRUE)
  }
  if (is.na(pkgfunlist[j, "fit"])) {
    return(rep(NA, 7))
  }
  if (echo > 1) {
    cat("____________________________________________________________________\n",
        j, pkgfunlist[j, "package"], "\n")
  } else if (echo > 0){
    cat("\n", j, pkgfunlist[j, "package"], "\n")
  }
  # Load package
  # do.call("library", list(pkgfunlist[j, "package"]))
  # Extract options
  other <- pkgotherpar[[pkgfunlist[j, "package"]]]

    fitargs <- list(obs)
    names(fitargs) <- pkgfunlist[j, "argdata"]
    res <- try(do.call(what =
                          getfun(paste0(list(pkgfunlist[j, "package"]),
                                        "::",
                                        pkgfunlist[j, "fit"])),
                       args = c(fitargs, other)),
               silent = TRUE)
   if(pkgfunlist[j,"package"] == "climextRemes") {
     # browser()
    if(res$info$failure) {
      res <- NULL
    }
  }

  if (inherits(res, "try-error") || is.null(res)) {
     res <- rep(NA, 7)
  } else {
    if(echo > 1) {
      cat("result before selection\n")
      print(res)
    }
    if(isS4(res)){
      res <- slot(res, pkgfunlist[j, "outputS4"])
    }
    if(!is.na(pkgfunlist[j, "outputS3"]))
      res <- res[[pkgfunlist[j, "outputS3"]]]
    if(!is.na(pkgfunlist[j, "outputS3par"]))
      res <- res[[pkgfunlist[j, "outputS3par"]]]
    if(echo > 1){
      cat("result after selection\n")
      print(res)
    }
      if(pkgfunlist[j,"package"] == "texmex") {
        # vector of coefficient with sigma = exp(phi), i.e.
        # optimization is performed on the log scale
        res <- c(res["mu: "], exp(res["phi: "]), res["xi: "])
        names(res)[1:3] <- c("mu", "sigma", "xi")
      } else if(!is.null(names(res))) {
        idxshape <- grep(pkgfunlist[j,"shape"], names(res), ignore.case = TRUE)
        idxscale <- grep(pkgfunlist[j,"scale"], names(res), ignore.case = TRUE)
        idxloc <- grep(pkgfunlist[j,"location"], names(res), ignore.case = TRUE)

        if(length(idxshape) > 0)
          names(res)[idxshape] <- pkgfunlist[j,"shape"]

        if(length(idxscale) > 0)
          names(res)[idxscale] <- pkgfunlist[j,"scale"]
        if(length(idxloc) > 0)
          names(res)[idxloc] <- pkgfunlist[j,"location"]
      }
    # if(pkgfunlist[j,"package"] == "climextRemes") {
    #    print(res)
    # }
      #reorder
      if(!is.null(names(res)) && all(pkgfunlist[j,c("location", "scale", "shape")] %in% names(res)))
        res <- res[c(pkgfunlist[j,"location"], pkgfunlist[j,"scale"], pkgfunlist[j,"shape"])]
      if(pkgfunlist[j, "package"] == "EnvStats"){
         res[3] <- -res[3] #shape parameter is flipped
      }
      # Compute log-likelihood
      ll <- try(mev::gev.ll(par = res[1:3], dat = obs),
                silent = TRUE)
      if (inherits(ll, "try-error")){
         ll <- NA
      }
      grad <- try(mev::gev.score(par = res[1:3],
                                 dat = obs),
                     silent = TRUE)
      if (inherits(grad, "try-error")){
         grad <- rep(NA, 3L)
      }
      res <- c(res, ll, grad)
      names(res) <- c("location","scale","shape","loglik",
                      "gradlocation","gradscale","gradshape")

  }
  if(is.list(res))
    res <- unlist(res)
  if(echo > 0)  {
    cat("MLE optimization done\n")
    print(res)
    cat("\n")
  }
  # mydetach(j, pkgfunlist)
  res
}

check_gev_fit_all <- function(obs,
                              pkgfunlist,
                              pkgotherpar,
                              ...){
  res <-
    sapply(
      1:NROW(pkgfunlist),
      check_gev_fit,
      obs = obs,
      pkgfunlist = pkgfunlist,
      pkgotherpar = pkgotherpar,
      ...
    )
  if (is.list(res)) {
    names(res) <- pkgfunlist[, "package"]

    if(all(sapply(res, length) == 7)) {
      res <- do.call(rbind, res)
    } else {
      print(res)
      stop("wrong result by check_gev_fit_all()")
    }
  } else {
    colnames(res) <- pkgfunlist[, "package"]
  }
  res
}

check_gev_varyingsize <- function(n = 10 ^ (3:4),
           nbrep = 2,
           dist,
           blocksize = 10,
           pkgfunlist,
           pkgotherpar,
           echo = FALSE,
           ...){
    nlen <- length(n)
    stopifnot(length(nbrep) == 1L)
    if(is.logical(echo)){
       echo <- as.integer(echo)
    }
    blocksize <- sort(as.integer(blocksize))
    stopifnot(length(blocksize) > 0)
    nb <- length(blocksize)
    res <- array(data = NA,
                 dim = c(nbrep, nb, nlen, 7, NROW(pkgfunlist)))
      dimnames(res) <-
        list(rep = paste0("rep", 1:nbrep),
             blocksize = paste0("b", blocksize),
             sampsize = paste0("n", n),
             param = c("location", "scale", "shape", "loglik",
               "gradlocation", "gradscale", "gradshape"),
             package = pkgfunlist[, "package"])

    for (nrep in seq_len(nbrep)) {
      if(echo > 0)
        cat("begin nrep", nrep, "\n")
      for (i in seq_len(nlen)) {
        if(echo > 0) {
          cat("*", i, "n_i", n[i], "*\n")
        }
         for (j in seq_along(blocksize)) {
             nbsimu <- blocksize[j] * n[i]
             obs <- matrix(do.call(paste0("r", dist),
                                   c(list(nbsimu), list(...))),
                           nrow = blocksize[j],
                           ncol = n[i])
             obs <- apply(obs, 2, max, na.rm = TRUE)

          temp <- check_gev_fit_all(obs,
                                    pkgfunlist,
                                    pkgotherpar,
                                    echo = echo)

             if(is.array(temp)) {
               # Return parameter estimates
               res[nrep, j, i, ,] <- temp
               thresh.mat <- NULL
             } else {
               print(class(temp))
               print(mode(temp))
               stop("wrong result")
             }
         }
         if(echo > 0) {
           cat("\n___end nrep", nrep, "___\n")
         }
      }
    }
   res[,,,,-which(apply(res[,1,1,1, ], 2, function(x) {
      sum(is.na(x))
   }) == nrep)]
}
