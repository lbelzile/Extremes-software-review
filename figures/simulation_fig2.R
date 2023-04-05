setwd(here::here())
source("../vignettes/1-univariate/utility-fun.R")

#### GPD list ####
gpdpkglist <- cbind(
  package =c("evd", "fExtremes", "evir", "extraDistr", "QRM",
             "ReIns", "TLMoments", "revdbayes", "Renext", "extRemes",
             "ismev", "lmomco", "SpatialExtremes", "texmex", "POT",
             "mev", "eva", "ercv", "lmom", "tea", "qrmtools") ,
  fun = c("gpd", "gpd", "gpd", "gpd", "GPD",
          "gpd", "gpd", "gp", "GPD", "evd",
          NA, "gpa", "gpd", "gpd", "gpd",
          NA,"gpd", NA, "gpa", "gpd", "GPD"),
  location = c("loc", "mu", "mu", "mu", NA,
               "mu", "loc", "loc", "loc", "loc",
               NA, "xi", "loc", "u", "loc",
               NA, "loc", NA, "xi", "loc","loc"),
  scale = c("scale", "beta", "beta", "sigma", "beta",
            "sigma", "scale", "scale", "scale", "scale",
            NA, "alpha", "scale", "sigma", "scale",
            "scale","scale", NA, "alpha", "scale","scale"),
  shape = c("shape", "xi", "xi", "xi", "xi",
            "gamma", "shape", "shape", "shape", "shape",
            NA, "kappa", "shape", "xi", "shape",
            "shape","shape", NA, "k", "shape", "shape"),
  fit = c("fpot", "gpdFit", "gpd", NA, "fit.GPD",
          "GPDfit", NA, NA, "fGPD", "fevd",
          "gpd.fit", NA, "gpdmle", "evm", "fitgpd",
          "fit.gpd","gpdFit", "fitpot", NA, "gpdFit","fit_GPD_MLE"),
  argdata = c("x", "x", "data", NA, "data",
              "data", NA, NA, "x", "x",
              "xdat", NA, "x", "y", "data",
              "xdat", "data", "data", NA, "data","x"),
  argthres = c("threshold","u", "threshold", NA, "threshold",
               NA, NA, NA, NA, "threshold",
               "threshold", NA, "threshold", "th", "threshold",
               "threshold","threshold", "threshold", NA, "threshold", NA),
  outputS4 = c(NA, "fit", NA, NA, NA,
               NA, NA, NA, NA, NA,
               NA, NA, NA, NA, NA,
               NA, NA, NA, NA, NA, NA),
  outputS3 = c("estimate", "fit", "par.ests", NA, "par.ests",
               NA, NA, NA, "estimate", "results",
               "mle", NA, NA, "coefficients", "fitted.values",
               "estimate","par.ests", "coeff", NA, "par.ests","par"),
  outputS3par = c(NA, "par", NA, NA, NA,
                  NA, NA, NA, NA, "par",
                  NA, NA, NA, NA, NA,
                  NA, NA, NA, NA, NA, NA)
)

rownames(gpdpkglist) <- gpdpkglist[,"package"]


mkfun <- function(x1, x2)
  ifelse(is.na(x2), NA, paste0(x1, x2))
gpdpkglist <- cbind(gpdpkglist, "dfun" = mkfun("d", gpdpkglist[,"fun"]))
gpdpkglist <- cbind(gpdpkglist, "pfun" = mkfun("p", gpdpkglist[,"fun"]))
#exception
gpdpkglist["lmomco", "pfun"] <- "cdfgpa"
gpdpkglist["lmomco", "dfun"] <- "pdfgpa"
gpdpkglist["lmom", "pfun"] <- "cdfgpa"
gpdpkglist["lmom", "dfun"] <- NA

#order
gpdpkglist <- gpdpkglist[order(gpdpkglist[,"package"]),]


gpdothpar.dist <- list(
  extRemes=list("type"="GP"))

gpdothpar.fit <- list(
  climextRemes = list("getParams"=TRUE),
  evd= list("std.err" = FALSE),
  extRemes=list("type"="GP"),
  ismev=list("show" = FALSE),
  QRM = list("verbose" = FALSE),
  POT = list("std.err.type"="none", "warn.inf"=FALSE),
  texmex=list("cov" = "numeric", "verbose"=FALSE),
  qrmtools = list("estimate.cov" = FALSE))

## ----install_packages----------------------------------
if(any(!gpdpkglist[,"package"] %in% installed.packages()[,"Package"]))
{
  pkg2install <- gpdpkglist[ !gpdpkglist[,"package"] %in% installed.packages()[,"Package"] ,"package"]
  cat(pkg2install, "\n")
  stop("missing packages")
}

## ---- GPD fit heavy tail -------------------------------------------
set.seed(123)
nrep <- 1000

time.gamma <- system.time(
  res1gamma <- check_gpd_gev_varyingsize(
  n = c(seq(from = 200L, to = 900L, by = 100L),
        seq(from = 1000L, to = 9000L, by = 1000L)),
  nbrep = 1000L,
  dist = "gamma",
  thres.qu = qgamma(p = 0.95, shape = 3, scale = 2),
  shape = 3,
  scale = 2,
  pkgfunlist = gpdpkglist,
  pkgotherpar = gpdothpar.fit,
  type = "gpd"))
save(res1gamma, file="Simulation_gamma.RData")
