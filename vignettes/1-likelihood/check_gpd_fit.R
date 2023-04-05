setwd(here::here())
source("utility-fun.R")
if (!dir.exists("outputs")) {
  dir.create("outputs")
}
cat("****** GPD test ******\n\n")

#### GPD list ####
gpdpkglist <- cbind(
  package = c(
    "evd",
    "fExtremes",
    "evir",
    "extraDistr",
    "QRM",
    "ReIns",
    "TLMoments",
    "revdbayes",
    "Renext",
    "extRemes",
    "ismev",
    "lmomco",
    "SpatialExtremes",
    "texmex",
    "POT",
    "mev",
    "eva",
    "ercv",
    "lmom",
    "tea",
    "qrmtools"
  ) ,
  fun = c(
    "gpd",
    "gpd",
    "gpd",
    "gpd",
    "GPD",
    "gpd",
    "gpd",
    "gp",
    "GPD",
    "evd",
    NA,
    "gpa",
    "gpd",
    "gpd",
    "gpd",
    NA,
    "gpd",
    NA,
    "gpa",
    "gpd",
    "GPD"
  ),
  location = c(
    "loc",
    "mu",
    "mu",
    "mu",
    NA,
    "mu",
    "loc",
    "loc",
    "loc",
    "loc",
    NA,
    "xi",
    "loc",
    "u",
    "loc",
    NA,
    "loc",
    NA,
    "xi",
    "loc",
    "loc"
  ),
  scale = c(
    "scale",
    "beta",
    "beta",
    "sigma",
    "beta",
    "sigma",
    "scale",
    "scale",
    "scale",
    "scale",
    NA,
    "alpha",
    "scale",
    "sigma",
    "scale",
    "scale",
    "scale",
    NA,
    "alpha",
    "scale",
    "scale"
  ),
  shape = c(
    "shape",
    "xi",
    "xi",
    "xi",
    "xi",
    "gamma",
    "shape",
    "shape",
    "shape",
    "shape",
    NA,
    "kappa",
    "shape",
    "xi",
    "shape",
    "shape",
    "shape",
    NA,
    "k",
    "shape",
    "shape"
  ),
  fit = c(
    "fpot",
    "gpdFit",
    "gpd",
    NA,
    "fit.GPD",
    "GPDfit",
    NA,
    NA,
    "fGPD",
    "fevd",
    "gpd.fit",
    NA,
    "gpdmle",
    "evm",
    "fitgpd",
    "fit.gpd",
    "gpdFit",
    "fitpot",
    NA,
    "gpdFit",
    "fit_GPD_MLE"
  ),
  argdata = c(
    "x",
    "x",
    "data",
    NA,
    "data",
    "data",
    NA,
    NA,
    "x",
    "x",
    "xdat",
    NA,
    "x",
    "y",
    "data",
    "xdat",
    "data",
    "data",
    NA,
    "data",
    "x"
  ),
  argthres = c(
    "threshold",
    "u",
    "threshold",
    NA,
    "threshold",
    NA,
    NA,
    NA,
    NA,
    "threshold",
    "threshold",
    NA,
    "threshold",
    "th",
    "threshold",
    "threshold",
    "threshold",
    "threshold",
    NA,
    "threshold",
    NA
  ),
  outputS4 = c(
    NA,
    "fit",
    NA,
    NA,
    NA,
    NA,
    NA,
    NA,
    NA,
    NA,
    NA,
    NA,
    NA,
    NA,
    NA,
    NA,
    NA,
    NA,
    NA,
    NA,
    NA
  ),
  outputS3 = c(
    "estimate",
    "fit",
    "par.ests",
    NA,
    "par.ests",
    NA,
    NA,
    NA,
    "estimate",
    "results",
    "mle",
    NA,
    NA,
    "coefficients",
    "fitted.values",
    "estimate",
    "par.ests",
    "coeff",
    NA,
    "par.ests",
    "par"
  ),
  outputS3par = c(
    NA,
    "par",
    NA,
    NA,
    NA,
    NA,
    NA,
    NA,
    NA,
    "par",
    NA,
    NA,
    NA,
    NA,
    NA,
    NA,
    NA,
    NA,
    NA,
    NA,
    NA
  )
)

rownames(gpdpkglist) <- gpdpkglist[, "package"]


mkfun <- function(x1, x2)
  ifelse(is.na(x2), NA, paste0(x1, x2))
gpdpkglist <-
  cbind(gpdpkglist, "dfun" = mkfun("d", gpdpkglist[, "fun"]))
gpdpkglist <-
  cbind(gpdpkglist, "pfun" = mkfun("p", gpdpkglist[, "fun"]))
#exception
gpdpkglist["lmomco", "pfun"] <- "cdfgpa"
gpdpkglist["lmomco", "dfun"] <- "pdfgpa"
gpdpkglist["lmom", "pfun"] <- "cdfgpa"
gpdpkglist["lmom", "dfun"] <- NA

#order
gpdpkglist <- gpdpkglist[order(gpdpkglist[, "package"]), ]

gpdpkglist <- cbind(gpdpkglist, "version" = sapply(gpdpkglist[,"package"], function(x){paste0(packageVersion(x), collapse = ".")}))

gpdothpar.dist <- list(extRemes = list("type" = "GP"))

gpdothpar.fit <- list(
  climextRemes = list("getParams" = TRUE),
  evd = list("std.err" = FALSE),
  extRemes = list("type" = "GP"),
  ismev = list("show" = FALSE),
  QRM = list("verbose" = FALSE),
  POT = list("std.err.type" = "none", "warn.inf" = FALSE),
  texmex = list("cov" = "numeric", "verbose" = FALSE),
  qrmtools = list("estimate.cov" = FALSE)
)

## ----install_packages----------------------------------
if (any(!gpdpkglist[, "package"] %in% installed.packages()[, "Package"]))
{
  pkg2install <-
    gpdpkglist[!gpdpkglist[, "package"] %in% installed.packages()[, "Package"] , "package"]
  cat(pkg2install, "\n")
  stop("missing packages")
}


rgpd <- function(n, scale, shape) {
  evd::rgpd(
    n = n,
    loc = 0,
    scale = scale,
    shape = shape
  )
}
n2 <- c(20L, 50L, 100L, 1000L)
fname <- "GPD-fit"

## ---- GPD fit heavy tail -------------------------------------------
set.seed(123)
nrep <- 1000


cat("\n******************************\n")
cat("***\t\t GPD fit 1 check\t\t***\n")


time.pos <- system.time(
  res_pos <-
    check_gpd_varyingsize(
      n = n2,
      nbrep = nrep,
      dist = "gpd",
      thres.qu = 0,
      shape = 0.5,
      scale = 1000,
      pkgfunlist = gpdpkglist,
      pkgotherpar = gpdothpar.fit,
      type = "gpd"
    )
)


save(res_pos, file = paste0(
  "outputs/",
  format(Sys.time(), "%Y_%m_%d_"),
  fname,
  "pos.RData"
))
print(time.pos)




## ---- GPD fit light tail -------------------------------------------
cat("\n******************************\n")
cat("***\t\t GPD fit 2 check\t\t***\n")
set.seed(123)
# qu <- qgamma(p = 0.95, shape = 3, scale = 2)
time.exp <- system.time(
  res_exp <- check_gpd_varyingsize(
    n = n2,
    nbrep = nrep,
    dist = "gpd",
    thres.qu = 0,
    scale = 1000,
    shape = 0,
    pkgfunlist = gpdpkglist,
    pkgotherpar = gpdothpar.fit,
    type = "gpd"
  )
)

save(res_exp, file = paste0("outputs/", format(Sys.time(), "%Y_%m_%d_"), fname, "exp.RData"))

## ---- GPD fit bounded tail -------------------------------------------


cat("\n******************************\n")
cat("***\t\t GPD fit 3 check\t\t***\n")

time.neg <- system.time(
  res_neg <- check_gpd_varyingsize(
    n = n2,
    nbrep = nrep,
    dist = "gpd",
    #thres.prob = 0.95,
    thres.qu = 0,
    shape = -0.5,
    scale = 1000,
    pkgfunlist = gpdpkglist,
    pkgotherpar = gpdothpar.fit,
    type = "gpd"
  )
)

save(res_neg, file = paste0("outputs/", format(Sys.time(), "%Y_%m_%d_"), fname, "neg.RData"))

write.csv(
  rbind(time.pos, time.exp, time.neg),
  paste0(
    "outputs/",
    format(Sys.time(), "%Y_%m_%d_"),
    "GPD-fitting-time.csv"
  )
)

cat("****** GPD finished ******\n\n")
