if (!"this.path" %in% installed.packages()[, "Package"]){
  install.packages("this.path")
}
setwd(this.path::here())
nrep <- 1000L
source("utility-fun-maxblock.R")
if(!dir.exists("outputs")){
  dir.create("outputs")
}
cat("****** GEV test ******\n\n")


#### GEV list ####

gevpkglist <- cbind(
  package =c("climextRemes", "evd", "evir","extraDistr", "extRemes",
             "fExtremes", "ismev", "lmomco", "QRM", "revdbayes",
             "SpatialExtremes", "texmex", "TLMoments","mev","qrmtools",
             "EnvStats","ExtremalDep") ,
  fun = c(NA, "gev", "gev", "gev", "evd",
          "gev", NA, "gev", "GEV", "gev",
          "gev", "gev", "gev", NA,"GEV", "gevd","GEV"),
  location = c("location", "loc", "mu", "mu", "loc",
               "mu", NA, "xi", "mu", "loc",
               "loc", "mu", "loc","loc","loc", "location","loc"),
  scale = c("scale", "scale", "sigma", "sigma", "scale",
            "beta", NA, "alpha", "sigma", "scale",
            "scale", "sigma", "scale", "scale","scale","scale","scale"),
  shape = c("shape", "shape", "xi", "xi", "shape", "xi", NA, "kappa",
            "xi", "shape", "shape", "xi", "shape", "shape","shape","shape","shape"),
  fit = c("fit_gev", "fgev", "gev", NA, "fevd",
          "gevFit", "gev.fit", NA, "fit.GEV", NA,
          "gevmle", "evm", NA, "fit.gev","fit_GEV_MLE","egevd","fGEV"),
  argdata = c("y", "x", "data", NA, "x",
              "x", "xdat", NA, "maxima", NA,
              "x", "y", NA, "xdat","x","x","data"),
  outputS4 = c(NA, NA, NA, NA, NA,
               "fit", NA, NA, NA, NA,
               NA, NA, NA, NA,NA, NA, NA),
  outputS3 = c("mle", "estimate", "par.ests", NA, "results",
               "par.ests", "mle", NA, "par.ests", NA,
               NA, "coefficients", NA, "param","par", "parameters", "est"),
  outputS3par = c(NA, NA, NA, NA, "par",
                  NA, NA, NA, NA, NA,
                  NA, NA, NA, NA, NA, NA))

gevothpar.dist <- list(
  extRemes=list("type"="GEV"))

gevothpar.fit <- list(
  climextRemes = list("getParams"=TRUE),
  evd= list(std.err = FALSE),
  extRemes=list("type"="GEV"),
  ismev=list("show" = FALSE),
  texmex=list("family"=texmex::gev, "verbose"=FALSE),
  qrmtools = list(estimate.cov = FALSE))

mkfun <- function(x1, x2)
  ifelse(is.na(x2), NA, paste0(x1, x2))

rownames(gevpkglist) <- gevpkglist[,"package"]

gevpkglist <- cbind(gevpkglist, "dfun" = mkfun("d", gevpkglist[,"fun"]))
gevpkglist <- cbind(gevpkglist, "pfun" = mkfun("p", gevpkglist[,"fun"]))

#exception
gevpkglist["lmomco", "pfun"] <- "cdfgev"
gevpkglist["lmomco", "dfun"] <- "pdfgev"
# Add package version
gevpkglist <- cbind(gevpkglist, "version" = sapply(gevpkglist[,"package"], function(x){paste0(packageVersion(x), collapse = ".")}))

# write.csv(gevpkglist, paste0("outputs/",format(Sys.time(), "%Y_%m_%d_%Hh%Mm%Ss_"), "gevpkglist-init.csv"))


## ----install_packages----------------------------------
if(any(!gevpkglist[,"package"] %in% installed.packages()[,"Package"])){
  pkg2install <- gevpkglist[ !gevpkglist[,"package"] %in% installed.packages()[,"Package"] ,"package"]
  cat(pkg2install, "\n")
  stop("missing packages")
}


cat("\n******************************\n")
cat("***\t\t GEV fit 1 check\t\t***\n")


n2 <- c(20L, 50L, 100L, 1000L)

fname <- "GEV-fit"
rlgamma <- function(n, shape, rate){
   exp(rgamma(n = n, shape = shape, rate = rate))
}
rgev <- evd::rgev
set.seed(202304)
time_pos <- system.time(
   res_pos <- check_gev_varyingsize(n = n2,
                                          blocksize = 1L,
                                          nbrep = nrep,
                                          dist = "gev",
                                          shape = 0.5,
                                          scale = 1000,
                                          loc = 0,
                                          pkgfunlist = gevpkglist,
                                          pkgotherpar = gevothpar.fit,
                                          echo = FALSE))


save(res_pos, file=paste0("outputs/",format(Sys.time(), "%Y_%m_%d_"), fname, "-pos", ".RData"))
print(time_pos)





## ---- GEV fit light tail -------------------------------------------

cat("\n******************************\n")
cat("***\t\t GEV fit 2 check\t\t***\n")

set.seed(202304)
time_zer <- system.time(
  res_zer <- check_gev_varyingsize(n = n2,
                                   blocksize = 1L,#c(10L, 25L, 50L),
                                     nbrep = nrep,
                                     dist = "gev",
                                     shape = 0,
                                     scale = 1000,
                                     loc = 0,
                                     pkgfunlist = gevpkglist,
                                     pkgotherpar = gevothpar.fit,
                                     echo = FALSE))

save(res_zer, file=paste0("outputs/",format(Sys.time(), "%Y_%m_%d_"), fname, "-zer", ".RData"))

print(time_zer)



## ---- GEV fit bounded tail -------------------------------------------

cat("\n******************************\n")
cat("***\t\t GEV fit 3 check\t\t***\n")

set.seed(202304)
time_neg <- system.time(
  res_neg <- check_gev_varyingsize(n = n2,
                                    blocksize = 1L,
                                    nbrep = nrep,
                                    dist = "gev",
                                    shape = -0.5,
                                    scale = 1000,
                                    pkgfunlist = gevpkglist,
                                    pkgotherpar = gevothpar.fit,
                                    echo = FALSE))

save(res_neg, file=paste0("outputs/",format(Sys.time(), "%Y_%m_%d_"), fname, "-neg", ".RData"))
print(time_neg)


#
# write.csv(rbind(time.loggamma, time.gamma, time.beta),
#           paste0("outputs/",format(Sys.time(), "%Y_%m_%d_"), "GEV-fitting-time.csv"))

cat("****** GEV finished ******\n\n")

