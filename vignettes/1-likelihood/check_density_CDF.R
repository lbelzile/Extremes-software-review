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
    NA
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
gpdpkglist <- gpdpkglist[order(gpdpkglist[, "package"]),]


gpdothpar.dist <- list(extRemes = list("type" = "GP"))

gpdothpar.fit <- list(
  climextRemes = list("getParams" = TRUE),
  evd = list("std.err" = FALSE),
  extRemes = list("type" = "GP"),
  ismev = list("show" = FALSE),
  QRM = list("verbose" = FALSE),
  POT = list("std.err.type" = "none", "warn.inf" = FALSE),
  texmex = list("cov" = "numeric", "verbose" = FALSE),
  qrmtools = list(estimate.cov = FALSE)
)



#
#
# write.csv(gpdpkglist, paste0("outputs/",format(Sys.time(), "%Y_%m_%d_"), "gpdpkglist-init.csv"))



## ----install_packages----------------------------------
if (any(!gpdpkglist[, "package"] %in% installed.packages()[, "Package"]))
{
  pkg2install <-
    gpdpkglist[!gpdpkglist[, "package"] %in% installed.packages()[, "Package"] , "package"]
  cat(pkg2install, "\n")
  stop("missing packages")
}


cat("\n******************************\n")
cat("***\t\t GPD density check\t\t***\n")



## ---- GPD compute density ----
xval <- seq(-2, 10, by = .1)
res_orig_pos_dens <-
  check_distrib_all(
    xval = xval,
    par0 = c(0, 1, 1),
    checkletter = "d",
    pkgfunlist = gpdpkglist,
    pkgotherpar = gpdothpar.dist
  )
res_orig_zer_dens <-
  check_distrib_all(
    xval = xval,
    par0 = c(0, 1, 0),
    checkletter = "d",
    pkgfunlist = gpdpkglist,
    pkgotherpar = gpdothpar.dist
  )
res_orig_neg_dens <-
  check_distrib_all(
    xval = xval,
    par0 = c(0, 1, -0.5),
    checkletter = "d",
    pkgfunlist = gpdpkglist,
    pkgotherpar = gpdothpar.dist
  )

shift <- 2
scale <- 2
shape_neg <- -0.5
res_shift_pos_dens <-
  check_distrib_all(
    xval = xval,
    par0 = c(shift, scale, 1),
    checkletter = "d",
    pkgfunlist = gpdpkglist,
    pkgotherpar = gpdothpar.dist
  )
res_shift_zer_dens <-
  check_distrib_all(
    xval = xval,
    par0 = c(shift, scale, 0),
    checkletter = "d",
    gpdpkglist,
    gpdothpar.dist
  )
res_shift_neg_dens <-
  check_distrib_all(
    xval = xval,
    par0 = c(shift, scale, -0.5),
    checkletter = "d",
    pkgfunlist = gpdpkglist,
    pkgotherpar = gpdothpar.dist
  )



## ---- Comments on density implementation -------------------------------------

zero_below_loc <-
  apply(res_shift_pos_dens[res_shift_pos_dens[, "x"] < shift, -1], 2, function(x) {
    ifelse(all(is.na(x)), NA,
           isTRUE(all.equal(max(abs(
             x
           )), 0, check.attributes = FALSE)))
  })
at_loc <-
  sapply(res_shift_pos_dens[res_shift_pos_dens[, "x"] == shift, -1], function(x) {
    ifelse(is.na(x), NA, isTRUE(all.equal(x, 0.5, check.attributes = FALSE)))
  })
zero_beyond_endpt <-
  apply(res_shift_neg_dens[res_shift_neg_dens[, "x"] > shift - scale / shape_neg, -1], 2, function(x) {
    ifelse(all(is.na(x)), NA, isTRUE(all.equal(max(abs(
      x
    )), 0, check.attributes = FALSE)))
  })
code_dens_valid <- ifelse(!zero_below_loc & !at_loc,
                          1,
                          ifelse(!zero_below_loc, 2,
                                 ifelse(!at_loc, 3,
                                        ifelse(is.na(
                                          at_loc
                                        ), 5, 4))))
gpdpkglist <- as.data.frame(gpdpkglist)
gpdpkglist$"density comment" <-
  c("incorrect for x < = loc",
    "incorrect for x < loc",
    "incorrect for x = loc",
    "correct",
    NA)[code_dens_valid]


gpd_pdf <- list(
  gpdpkglist,
  "orig_pos" = res_orig_pos_dens,
  "orig_zero" = res_orig_zer_dens,
  "orig_neg" = res_orig_neg_dens,
  "shift_pos" = res_shift_pos_dens,
  "shift_zero" = res_shift_zer_dens,
  "shift_neg" = res_shift_neg_dens,
  "zero_below_loc" = zero_below_loc,
  "at_loc" = at_loc,
  "zero_beyond_endpt" = zero_beyond_endpt
)

save(
  gpd_pdf,
  file = paste0(
    "outputs/",
    format(Sys.time(), "%Y_%m_%d_"),
    "GPD-density-function.RData"
  ),
  version = 2
)
cat("end GPD density function test\n")


## ---- GPD compute distribution function ----


cat("\n******************************\n")
cat("***\t\t GPD cumul distrib check\t\t***\n")

xval <- seq(0, 10, by = .1)
res_orig_pos_cdf <-
  check_distrib_all(
    xval = xval,
    par0 = c(0, 1, 1),
    checkletter = "p",
    pkgfunlist = gpdpkglist,
    pkgotherpar = gpdothpar.dist
  )
res_orig_zer_cdf <-
  check_distrib_all(
    xval = xval,
    par0 = c(0, 1, 0),
    checkletter = "p",
    pkgfunlist = gpdpkglist,
    pkgotherpar = gpdothpar.dist
  )
res_orig_neg_cdf <-
  check_distrib_all(
    xval = xval,
    par0 = c(0, 1, -0.5),
    checkletter = "p",
    pkgfunlist = gpdpkglist,
    pkgotherpar = gpdothpar.dist
  )



xval <- seq(-1, 10, by = .1)
res_shift_pos_cdf <-
  check_distrib_all(
    xval = xval,
    par0 = c(0, 2, 1),
    checkletter = "p",
    pkgfunlist = gpdpkglist,
    pkgotherpar = gpdothpar.dist
  )
res_shift_zer_cdf <-
  check_distrib_all(
    xval = xval,
    par0 = c(0, 2, 0),
    checkletter = "p",
    gpdpkglist,
    gpdothpar.dist
  )
res_shift_neg_cdf <-
  check_distrib_all(
    xval = xval,
    par0 = c(0, 2, -0.5),
    checkletter = "p",
    pkgfunlist = gpdpkglist,
    pkgotherpar = gpdothpar.dist
  )



## ----Comments on distribution function implementation ------------------

zero_below_loc <-
  apply(res_shift_neg_cdf[res_shift_neg_cdf[, "x"] < 0, -1],
        2,
        function(x) {
          ifelse(all(is.na(x)), NA, isTRUE(all.equal(max(abs(
            x
          )), 0, check.attributes = FALSE)))
        })
zero_at_loc <-
  sapply(res_shift_neg_cdf[res_shift_neg_cdf[, "x"] == 0, -1], function(x) {
    ifelse(all(is.na(x)), NA, isTRUE(all.equal(max(abs(
      x
    )), 0, check.attributes = FALSE)))
  })
one_at_endpoint <-
  sapply(res_shift_neg_cdf[res_shift_neg_cdf[, "x"] == 8, -1],  function(x) {
    ifelse(is.na(x), NA, isTRUE(all.equal(x, 1, check.attributes = FALSE)))
  })
one_beyond_endpoint <-
  apply(res_shift_neg_cdf[res_shift_neg_cdf[, "x"] > 6, -1], 2, function(x) {
    ifelse(all(is.na(x)), NA, isTRUE(all.equal(min(x), 1, check.attributes = FALSE)))
  })

code_cdf_valid <- ifelse(!zero_below_loc | !one_beyond_endpoint,
                         1,
                         ifelse(!zero_below_loc, 2,
                                ifelse(
                                  !one_beyond_endpoint, 6,
                                  ifelse(!zero_at_loc, 3,
                                         ifelse(is.na(zero_at_loc), 5, 4))
                                )))
gpdpkglist <- as.data.frame(gpdpkglist)
gpdpkglist$"distribution comment" =
  c(
    "incorrect outside support",
    "incorrect for x < loc",
    "incorrect for x = loc",
    "correct",
    NA,
    "incorrect beyond endpoint"
  )[code_cdf_valid]

gpd_cdf <-
  list(
    "orig_pos" = res_orig_pos_cdf,
    "orig_zero" = res_orig_zer_cdf,
    "orig_neg" = res_orig_neg_cdf,
    "shift_pos" = res_shift_pos_cdf,
    "shift_zero" = res_shift_zer_cdf,
    "shift_neg" = res_shift_neg_cdf,
    "zero_below_loc" = zero_below_loc,
    "zero_at_loc" = zero_at_loc,
    "one_at_endpoint" = one_at_endpoint,
    "one_beyond_endpoint" = one_beyond_endpoint
  )


save(
  gpd_cdf,
  file = paste0(
    "outputs/",
    format(Sys.time(), "%Y_%m_%d_"),
    "GPD-distrib-function.RData"
  ),
  version = 2
)

write.csv(gpdpkglist,
          paste0(
            "outputs/",
            format(Sys.time(), "%Y_%m_%d_"),
            "gpdpkglist-final.csv"
          ),
          row.names = FALSE)
cat("end GPD distrib function test\n")




source("utility-fun-maxblock.R")


#### GEV list ####

gevpkglist <- cbind(
  package = c(
    "climextRemes",
    "evd",
    "evir",
    "extraDistr",
    "extRemes",
    "fExtremes",
    "ismev",
    "lmomco",
    "QRM",
    "revdbayes",
    "SpatialExtremes",
    "texmex",
    "TLMoments",
    "qrmtools",
    "EnvStats",
    "ExtremalDep"
  ) ,
  fun = c(
    NA,
    "gev",
    "gev",
    "gev",
    "evd",
    "gev",
    NA,
    "gev",
    "GEV",
    "gev",
    "gev",
    "gev",
    "gev",
    "GEV",
    "gevd",
    "GEV"
  ),
  location = c(
    "location",
    "loc",
    "mu",
    "mu",
    "loc",
    "mu",
    NA,
    "xi",
    "mu",
    "loc",
    "loc",
    "mu",
    "loc",
    "loc",
    "location",
    "loc"
  ),
  scale = c(
    "scale",
    "scale",
    "sigma",
    "sigma",
    "scale",
    "beta",
    NA,
    "alpha",
    "sigma",
    "scale",
    "scale",
    "sigma",
    "scale",
    "scale",
    "scale",
    "scale"
  ),
  shape = c(
    "shape",
    "shape",
    "xi",
    "xi",
    "shape",
    "xi",
    NA,
    "kappa",
    "xi",
    "shape",
    "shape",
    "xi",
    "shape",
    "shape",
    "shape",
    "shape"
  ),
  fit = c(
    "fit_gev",
    "fgev",
    "gev",
    NA,
    "fevd",
    "gevFit",
    "gev.fit",
    NA,
    "fit.GEV",
    NA,
    NA,
    "evm",
    NA,
    "fit_GEV_MLE",
    "egevd",
    "fGEV"
  ),
  argdata = c(
    "y",
    "x",
    "data",
    NA,
    "x",
    "x",
    "xdat",
    NA,
    "maxima",
    NA,
    NA,
    "y",
    NA,
    "x",
    "x",
    "data"
  ),
  outputS4 = c(NA, NA, NA, NA, NA,
               "fit", NA, NA, NA, NA,
               NA, NA, NA, NA, NA, NA),
  outputS3 = c(
    "mle",
    "estimate",
    "par.ests",
    NA,
    "results",
    "par.ests",
    "mle",
    NA,
    "par.ests",
    NA,
    NA,
    "coefficients",
    NA,
    "par",
    "parameters",
    "est"
  ),
  outputS3par = c(NA, NA, NA, NA, "par",
                  NA, NA, NA, NA, NA,
                  NA, NA, NA, NA, NA, NA)
)

gevothpar.dist <- list(extRemes = list("type" = "GEV"))

gevothpar.fit <- list(
  climextRemes = list("getParams" = TRUE),
  evd = list(std.err = FALSE),
  extRemes = list("type" = "GEV"),
  ismev = list("show" = FALSE),
  texmex = list("family" = texmex::gev, "verbose" = FALSE),
  qrmtools = list(estimate.cov = FALSE)
)

mkfun <- function(x1, x2)
  ifelse(is.na(x2), NA, paste0(x1, x2))

rownames(gevpkglist) <- gevpkglist[, "package"]

gevpkglist <-
  cbind(gevpkglist, "dfun" = mkfun("d", gevpkglist[, "fun"]))
gevpkglist <-
  cbind(gevpkglist, "pfun" = mkfun("p", gevpkglist[, "fun"]))

#exception
gevpkglist["lmomco", "pfun"] <- "cdfgev"
gevpkglist["lmomco", "dfun"] <- "pdfgev"


# write.csv(gevpkglist, paste0("outputs/",format(Sys.time(), "%Y_%m_%d_"), "gevpkglist-init.csv"))


## ----install_packages----------------------------------
if (any(!gevpkglist[, "package"] %in% installed.packages()[, "Package"]))
{
  pkg2install <-
    gevpkglist[!gevpkglist[, "package"] %in% installed.packages()[, "Package"] , "package"]
  cat(pkg2install, "\n")
  stop("missing packages")
}


cat("\n******************************\n")
cat("***\t\t GEV density check\t\t***\n")


## ---- GEV compute density ----
xval <- seq(-5, 10, by = .1)
res_orig_pos_dens <-
  check_distrib_all(
    xval = xval,
    par0 = c(0, 1, 1),
    checkletter = "d",
    pkgfunlist = gevpkglist,
    gevothpar.dist
  )

res_orig_zer_dens <-
  check_distrib_all(
    xval = xval,
    par0 = c(0, 1, 0),
    checkletter = "d",
    pkgfunlist = gevpkglist,
    gevothpar.dist
  )

res_orig_almostzer_dens <-
  check_distrib_all(
    xval = xval,
    par0 = c(0, 1, 1e-8),
    checkletter = "d",
    pkgfunlist = gevpkglist,
    gevothpar.dist
  )

res_orig_neg_dens <-
  check_distrib_all(
    xval = xval,
    par0 = c(0, 1, -1 / 2),
    checkletter = "d",
    pkgfunlist = gevpkglist,
    gevothpar.dist
  )

xval <- seq(-5, 10, by = .1)

res_shift_pos_dens <-
  check_distrib_all(
    xval = xval,
    par0 = c(2, 2, 1),
    checkletter = "d",
    pkgfunlist = gevpkglist,
    gevothpar.dist
  )

res_shift_zer_dens <-
  check_distrib_all(
    xval = xval,
    par0 = c(2, 2, 0),
    checkletter = "d",
    pkgfunlist = gevpkglist,
    gevothpar.dist
  )

res_shift_neg_dens <-
  check_distrib_all(
    xval = xval,
    par0 = c(2, 2, -1 / 2),
    checkletter = "d",
    pkgfunlist = gevpkglist,
    gevothpar.dist
  )




## ---- Comments on density implementation -------------------------------------

zero_below_loc <-
  apply(res_shift_pos_dens[res_shift_pos_dens[, "x"] < 2 - 2 / 1, -1], 2, function(x) {
    ifelse(all(is.na(x)), NA, isTRUE(all.equal(max(abs(
      x
    )), 0, check.attributes = FALSE)))
  })
at_loc <-
  sapply(res_shift_pos_dens[res_shift_pos_dens[, "x"] == 2 - 2 / 1, -1], function(x) {
    ifelse(is.na(x), NA, isTRUE(!is.nan(x)))
  })
zero_beyond_endpt <-
  apply(res_shift_neg_dens[res_shift_neg_dens[, "x"] > 2 - 2 / (-1 / 2), -1], 2, function(x) {
    ifelse(all(is.na(x)), NA, isTRUE(all.equal(max(abs(
      x
    )), 0, check.attributes = FALSE)))
  })
code_dens_valid <- ifelse(!zero_below_loc & !at_loc,
                          1,
                          ifelse(!zero_below_loc, 2,
                                 ifelse(!at_loc, 3,
                                        ifelse(is.na(
                                          at_loc
                                        ), 5, 4))))
gevpkglist <- as.data.frame(gevpkglist)
gevpkglist$"density comment" <-
  c("incorrect for x < = loc",
    "incorrect for x < loc",
    "incorrect for x = loc",
    "correct",
    NA)[code_dens_valid]


gev_pdf <-
  list(
    "orig_pos" = res_orig_pos_dens,
    "orig_zero" = res_orig_zer_dens,
    "orig_neg" = res_orig_neg_dens,
    "shift_pos" = res_shift_pos_dens,
    "shift_zero" = res_shift_zer_dens,
    "shift_neg" = res_shift_neg_dens,
    "zero_below_loc" = zero_below_loc,
    "at_loc" = at_loc,
    "zero_beyond_endpt" = zero_beyond_endpt
  )


save(
  gev_pdf,
  file = paste0(
    "outputs/",
    format(Sys.time(), "%Y_%m_%d_"),
    "GEV-density-function.RData"
  ),
  version = 2
)
cat("end GEV density function test\n")


## ---- GEV compute distribution function ----

cat("\n******************************\n")
cat("***\t\t GEV cumul distrib check\t\t***\n")


xval <- seq(-5, 10, by = .1)
res_orig_pos_cdf <-
  check_distrib_all(
    xval = xval,
    par0 = c(0, 1, 1),
    checkletter = "p",
    pkgfunlist = gevpkglist,
    gevothpar.dist
  )

res_orig_zer_cdf <-
  check_distrib_all(
    xval = xval,
    par0 = c(0, 1, 0),
    checkletter = "p",
    pkgfunlist = gevpkglist,
    gevothpar.dist
  )

res_orig_almostzer_cdf <-
  check_distrib_all(
    xval = xval,
    par0 = c(0, 1, 0),
    checkletter = "p",
    pkgfunlist = gevpkglist,
    gevothpar.dist
  )


res_orig_neg_cdf <-
  check_distrib_all(
    xval = xval,
    par0 = c(0, 1, -0.5),
    checkletter = "p",
    pkgfunlist = gevpkglist,
    gevothpar.dist
  )

xval <- seq(-5, 10, by = .1)

res_shift_pos_cdf <-
  check_distrib_all(
    xval = xval,
    par0 = c(2, 2, 1),
    checkletter = "p",
    pkgfunlist = gevpkglist,
    gevothpar.dist
  )

res_shift_zer_cdf <-
  check_distrib_all(
    xval = xval,
    par0 = c(2, 2, 0),
    checkletter = "p",
    pkgfunlist = gevpkglist,
    gevothpar.dist
  )

res_shift_almostzer_cdf <-
  check_distrib_all(
    xval = xval,
    par0 = c(2, 2, 1e-8),
    checkletter = "p",
    pkgfunlist = gevpkglist,
    gevothpar.dist
  )

res_shift_neg_cdf <-
  check_distrib_all(
    xval = xval,
    par0 = c(2, 2, -0.5),
    checkletter = "p",
    pkgfunlist = gevpkglist,
    gevothpar.dist
  )
#
# res_nearzero_cdf <-
#   check_distrib_all(
#     xval = xval,
#     par0 = c(2, 2, 1e-10),
#     checkletter = "d",
#     pkgfunlist = gevpkglist,
#     gevothpar.dist
#   )
# matplot(res_nearzero_cdf[,1], res_nearzero_cdf[,-1], type = "l")



## ---- Comments on distrib fun implementation -------------------------------------
#
zero_below_loc <-
  apply(res_shift_pos_cdf[res_shift_pos_cdf[, "x"] < 2 - 2 / 1, -1], 2, function(x) {
    ifelse(all(is.na(x)), NA, isTRUE(all.equal(max(abs(
      x
    )), 0, check.attributes = FALSE)))
  })
at_loc <-
  sapply(res_shift_pos_cdf[res_shift_pos_cdf[, "x"] == 2 - 2 / 1, -1], function(x) {
    ifelse(is.na(x), NA, isTRUE(!is.nan(x)))
  })
zero_beyond_endpt <-
  apply(res_shift_neg_cdf[res_shift_neg_cdf[, "x"] > 2 - 2 / (-1 / 2), -1], 2, function(x) {
    ifelse(all(is.na(x)), NA, isTRUE(all.equal(max(abs(
      x
    )), 0, check.attributes = FALSE)))
  })
code_cdf_valid <- ifelse(!zero_below_loc & !at_loc,
                         1,
                         ifelse(!zero_below_loc, 2,
                                ifelse(!at_loc, 3,
                                       ifelse(is.na(
                                         at_loc
                                       ), 5, 4))))
gevpkglist <- as.data.frame(gevpkglist)
gevpkglist$"distribution comment" <-
  c("incorrect for x < = loc",
    "incorrect for x < loc",
    "incorrect for x = loc",
    "correct",
    NA)[code_dens_valid]

gev_cdf <-
  list(
    "orig_pos" = res_orig_pos_cdf,
    "orig_zero" = res_orig_zer_cdf,
    "orig_neg" = res_orig_neg_cdf,
    "shift_pos" = res_shift_pos_cdf,
    "shift_zero" = res_shift_zer_cdf,
    "shift_neg" = res_shift_neg_cdf,
    "zero_below_loc" = zero_below_loc,
    "at_loc" = at_loc,
    "zero_beyond_endpt" = zero_beyond_endpt
  )

save(
  gev_cdf,
  file = paste0(
    "outputs/",
    format(Sys.time(), "%Y_%m_%d_"),
    "GEV-distrib-function.RData"
  ),
  version = 2
)

write.csv(
  gevpkglist,
  file = paste0(
    "outputs/",
    format(Sys.time(), "%Y_%m_%d_"),
    "gevpkglist-final.csv"
  ),
  row.names = FALSE
)
cat("end GEV distrib function test\n")
