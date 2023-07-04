# https://stackoverflow.com/questions/38983179/do-call-a-function-in-r-without-loading-the-package
getfun <- function(x) {
  if(length(grep("::", x))>0) {
    parts<-strsplit(x, "::")[[1]]
    getExportedValue(parts[1], parts[2])
  } else {
    x
  }
}


## ----internal_functions, echo = FALSE------------------
mydetach <- function(j, pkglist)
{
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
myfulldetach <- function()
{
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
myfullunload <- function()
{
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
check_distrib <-
  function(j,
           xval,
           par0,
           checkletter = "d",
           echo = FALSE,
           pkgfunlist,
           pkgotherpar=NULL)
  {
    if (is.na(pkgfunlist[j, "fun"]))
      return(rep(NA, length(xval)))
    if (is.na(pkgfunlist[j, "dfun"]) && checkletter == "d")
      return(rep(NA, length(xval)))
    if (is.na(pkgfunlist[j, "pfun"]) && checkletter == "p")
      return(rep(NA, length(xval)))

    cat("\n******", j, "******\n")

    # do.call("library", list(pkgfunlist[j, "package"]))

    if (!is.na(pkgfunlist[j, "location"]))
    {
      names(par0) <- c(pkgfunlist[j, "location"],
                       pkgfunlist[j, "scale"],
                       pkgfunlist[j, "shape"])
    } else
    {
      xval <- xval - par0[1]
      par0 <- par0[2:3]
      names(par0) <- c(pkgfunlist[j, "scale"],
                       pkgfunlist[j, "shape"])
    }
    other <- pkgotherpar[[pkgfunlist[j, "package"]]]
    if(pkgfunlist[j, "package"] == "EnvStats"){
       par0 <- par0[c(pkgfunlist[j, "location"],
                      pkgfunlist[j, "scale"],
                      pkgfunlist[j, "shape"])]
       par0[3] <- -par0[3] # parameterized in terms of negative shape
    }
    if (pkgfunlist[j, "package"] == "lmomco")
    {
      #reorder
      #par0[pkgfunlist[j, "shape"]] <- -par0[pkgfunlist[j, "shape"]]
      par0 <- par0[c(pkgfunlist[j, "location"],
                     pkgfunlist[j, "scale"],
                     pkgfunlist[j, "shape"])]
      #convert
      par0 <-
        list(lmomco::vec2par(
          vec = par0,
          type = pkgfunlist[j, "fun"],
          paracheck = FALSE
        ))
      names(par0) <- "para"
    } else if (pkgfunlist[j, "package"] == "lmom")
    {
      par0 <- list(as.numeric(par0))
      names(par0) <- "para"
    }
    if (echo)
    {
      cat(j, pkgfunlist[j, "package"], "\n")
    }
    if (checkletter == "d") {
      res <-
        do.call(getfun(paste0(list(pkgfunlist[j, "package"]),
                              "::",
                              pkgfunlist[j, "dfun"])),
                c(list(xval), as.list(par0), other))
    } else if (checkletter == "p") {
      res <-
        do.call(getfun(paste0(list(pkgfunlist[j, "package"]),
                              "::",
                              pkgfunlist[j, "pfun"])),
                c(list(xval), as.list(par0), other))
    } else{
      stop("wrong function")
    }
    # myfulldetach()
    cat("\n")
    #myfullunload() #unneeded
    as.numeric(res)
  }
check_distrib_all <- function(xval, par0, checkletter,
                              pkgfunlist, pkgotherpar)
{
  res <-
    cbind(
      xval,
      sapply(
        1:NROW(pkgfunlist),
        check_distrib,
        xval = xval,
        par0 = par0,
        checkletter = checkletter,
        pkgfunlist = pkgfunlist,
        pkgotherpar = pkgotherpar
      )
    )
  colnames(res) <- c("x", pkgfunlist[, "package"])
  res
}
check_gpd_fit <- function(j, obs, thres, echo = FALSE,
                          pkgfunlist, pkgotherpar, type="gpd")
{
  if(is.logical(echo))
    echo <- 1*echo
  type <- match.arg(type, c("gpd", "gev"))
  if (is.na(pkgfunlist[j, "fit"]))
  {
    if(type == "gpd")
      return(rep(NA, 2))
    else
      return(rep(NA, 3))
  }
  if (echo > 1){
    cat("____________________________________________________________________\n",
        j, pkgfunlist[j, "package"], "\n")
  } else if (echo > 0){
    cat("\n", j, pkgfunlist[j, "package"], "\n")
  }
  # do.call("library", list(pkgfunlist[j, "package"]))

  other <- pkgotherpar[[pkgfunlist[j, "package"]]]

  if(type == "gpd")
  {
    if (!is.na(pkgfunlist[j, "argthres"]))
    {
      fitargs <- list(obs, thres)
      names(fitargs) <- pkgfunlist[j, c("argdata", "argthres")]
      res <- try(do.call(what = getfun(paste0(list(pkgfunlist[j, "package"]),
                               "::",pkgfunlist[j, "fit"])),
                         args = c(fitargs, other)), silent = TRUE)
    } else
    {
      #compute excesses if no threshold argument
      fitargs <- list(obs[obs > thres] - thres)
      names(fitargs) <- pkgfunlist[j, "argdata"]
      res <- try(do.call(what = getfun(paste0(list(pkgfunlist[j, "package"]),
                                              "::",pkgfunlist[j, "fit"])),
                         args = c(fitargs, other)),
                 silent = TRUE)
    }
  } else
  {
    fitargs <- list(obs)
    names(fitargs) <- pkgfunlist[j, "argdata"]
    res <- try(do.call(what = getfun(paste0(list(pkgfunlist[j, "package"]),
                                            "::",
                                            pkgfunlist[j, "fit"])),
                       args = c(fitargs, other)), silent = TRUE)
  }
  if(pkgfunlist[j,"package"] == "climextRemes")
  {
    if(res$info$failure)
    {
      res <- NULL
    }
  }

  if (inherits(res, "try-error") || is.null(res))
  {
    if(type == "gpd")
    {
      res <- rep(NA, 2)
    } else
      res <- rep(NA, 3)
  } else
  {
    if(echo > 1)
    {
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
    if(echo > 1)
    {
      cat("result after selection\n")
      print(res)
    }
    if(type == "gpd")
    {
      #check names => reformat names if needed
      if(pkgfunlist[j,"package"] == "texmex")
      {
        # vector of coefficient with sigma = exp(phi), i.e.
        # optimization is performed on the log scale
        res <- c(res["xi: "], exp(res["phi: "]))
        names(res) <- c("shape", "scale")
      } else if(pkgfunlist[j,"package"] == "ercv")
      {
        res <- res[c("evi", "psi")]
        names(res) <- c("shape", "scale")
      } else if(pkgfunlist[j,"package"] == "ismev")
      {
        res <- c("shape" = res[2], "scale" = res[1])
      } else
      {
        if(!is.null(names(res)))
        {
          idxshape <- grep(pkgfunlist[j,"shape"], names(res), ignore.case = TRUE)
          idxscale <- grep(pkgfunlist[j,"scale"], names(res), ignore.case = TRUE)

          if(length(idxshape) > 0)
          {
            names(res)[idxshape] <- pkgfunlist[j,"shape"]
            if(length(idxscale) > 0)
              names(res)[idxscale] <- pkgfunlist[j,"scale"]
            else
              names(res)[-idxshape] <- pkgfunlist[j,"scale"]
          }
        }
        #reorder
        if(!is.null(names(res)))
          res <- res[c(pkgfunlist[j,"shape"], pkgfunlist[j,"scale"])]

      }
    } else #type == "gev"
    {
      if(pkgfunlist[j,"package"] == "texmex")
      {
        # vector of coefficient with sigma = exp(phi), i.e.
        # optimization is performed on the log scale
        res <- c(res["mu: "], exp(res["phi: "]), res["xi: "])
        names(res) <- c("mu", "sigma", "xi")
      } else if(!is.null(names(res)))
      {
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
      #reorder
      if(!is.null(names(res)) && all(pkgfunlist[j,c("location", "scale", "shape")] %in% names(res)))
        res <- res[c(pkgfunlist[j,"location"], pkgfunlist[j,"scale"], pkgfunlist[j,"shape"])]

    }
  }
  if(is.list(res))
    res <- unlist(res)
  if(echo > 0)
  {
    cat("MLE optimization done\n")
    print(res)
    cat("\n")
  }
  # mydetach(j, pkgfunlist)
  res
}

check_gpd_fit_all <- function(obs, thres=NULL, pkgfunlist,
                              pkgotherpar, type, ...)
{
  res <-
    sapply(
      1:NROW(pkgfunlist),
      check_gpd_fit,
      obs = obs,
      thres = thres,
      pkgfunlist = pkgfunlist,
      pkgotherpar = pkgotherpar,
      type = type,
      ...
    )
  if (is.list(res))
  {
    names(res) <- pkgfunlist[, "package"]

    if(all(sapply(res, length) == 2) && type == "gpd")
      res <- do.call(rbind, res)
    else if(all(sapply(res, length) == 3) && type == "gev")
      res <- do.call(rbind, res)
    else
    {
      print(res)
      stop("wrong result by check_gpd_fit_all()")
    }
  } else
    colnames(res) <- pkgfunlist[, "package"]
  res
}

check_gpd_varyingsize <-function(n = 10 ^ (3:4),
           nbrep = 2,
           dist,
           thres.prob,
           thres.qu = NULL,
           pkgfunlist,
           pkgotherpar,
           echo = FALSE,
           type,
           ...)
  {
    nlen <- length(n)
    if(is.logical(echo))
      echo <- 1*echo
    if(type == "gpd")
    {
      res <- array(data = NA, dim = c(nbrep, nlen, 5, NROW(pkgfunlist)))
      if(is.null(thres.qu)){
        thresh.mat <- array(dim = c(nbrep, nlen))
      } else{
        thresh.mat <- NA
      }
      dimnames(res) <-
        list(paste0("rep", 1:nbrep),
             paste0("n", n),
             c("shape", "scale", "loglik", "gradshape", "gradscale"),
             pkgfunlist[, "package"])
    } else if(type == "gev")
    {
      res <- array(data = NA, dim = c(nbrep, nlen, 7, NROW(pkgfunlist)))
      dimnames(res) <-
        list(paste0("rep", 1:nbrep),
             paste0("n", n),
             c("location", "scale", "shape", "loglik",
               "gradlocation", "gradscale", "gradshape"),
             pkgfunlist[, "package"])
    } else
      stop("wrong type")

    for (nrep in 1:nbrep)
    {
      if(echo > 0)
        cat("begin nrep", nrep, "\n")
      for (i in 1:nlen)
      {
        if(echo > 0)
        {
          cat("*", i, "n_i", n[i], "*\n")
        }
        obs <- do.call(paste0("r", dist), c(list(n[i]), list(...)))
        if(!is.null(thres.qu))
        {
          thres <- thres.qu[1]
        } else{
          thres <- quantile(obs, probs = thres.prob)
          thresh.mat[nrep,i] <- thres
        }

        if(type == "gpd")
        {
          temp <- check_gpd_fit_all(obs, thres, pkgfunlist,
                                    pkgotherpar, type, echo=echo)
          if(is.array(temp))
          {
            # Return parameter estimates (scale and shape)
            res[nrep, i, c("shape", "scale"),] <- temp
          } else
            stop("wrong result")

          # Log-likelihood using a correct implementation of the GP density
          exc <- obs[obs>thres] - thres
          for(j in seq_len(NROW(pkgfunlist)))
          {
            if(!is.na(res[nrep, i, 1,j]))
            {
              # Compute log-likelihood
              ll <- try(mev::gpd.ll(par = res[nrep, i, 2:1,j], dat = exc),
                        silent = TRUE)
              if (!inherits(ll, "try-error")){
                res[nrep, i, "loglik", j] <- ll
                grad <- try(mev::gpd.score(par = res[nrep, i, 2:1,j], dat = exc),
                            silent = TRUE)
                if (!inherits(grad, "try-error")){
                  res[nrep, i, c("gradshape", "gradscale"), j] <- rev(grad)
                }
              }

            }
          }
        } else if(type == "gev")
        {
          temp <- check_gpd_fit_all(obs, thres, pkgfunlist,
                                    pkgotherpar, type, echo=echo)

          if(is.array(temp))
          {
            # Return parameter estimates
          res[nrep, i, c("location", "scale", "shape"),] <- temp
          thresh.mat <- NULL
          } else
          {
            print(class(temp))
            print(mode(temp))
            stop("wrong result")
          }

          for(j in seq_len(NROW(pkgfunlist)))
          {
            if(!is.na(res[nrep, i, 1,j]))
            {
              # Compute log-likelihood
              ll <- try(mev::gev.ll(par = res[nrep, i, 2:1,j], dat = exc),
                        silent = TRUE)
              if (!inherits(ll, "try-error")){
                res[nrep, i, "loglik", j] <- ll
                grad <- try(mev::gev.score(par = res[nrep, i, 2:1,j], dat = exc),
                            silent = TRUE)
                if (!inherits(grad, "try-error")){
                  res[nrep, i, c("gradlocation", "gradscale", "gradshape"), j] <- grad
                }
              }

            }
          }
        }
      }
      if(echo > 0){
        cat("\n___end nrep", nrep, "___\n")
      }
    }
    res.avg <-
      apply(res[, , , pkgfunlist[!is.na(pkgfunlist[, "fit"]), "package"]], 2:4,
            mean, na.rm = TRUE)
    res.sd <-
      apply(res[, , , pkgfunlist[!is.na(pkgfunlist[, "fit"]), "package"]], 2:4,
            sd, na.rm = TRUE)
    res.noncvg <-
      apply(res[, , , pkgfunlist[!is.na(pkgfunlist[, "fit"]), "package"]], 2:4,
            function(x)
              sum(is.na(x)))

    list(
      fullres = res,
      avg = res.avg,
      sd = res.sd,
      noncvg = res.noncvg,
      thresh = thresh.mat
    )
  }
