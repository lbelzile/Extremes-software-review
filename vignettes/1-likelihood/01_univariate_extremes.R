#| label: setup
#| echo: false
#| eval: true
library(ggplot2)
library(patchwork)
if (!"this.path" %in% installed.packages()[, "Package"]){
  install.packages("this.path")
}
setwd(this.path::here())


#| label: tbl-gp1
#| tbl-cap: "Evaluation of  generalized Pareto model density and distribution functions."
#| eval: true
#| echo: false
#| message: false
#| warning: false
#| cache: true
options(knitr.kable.NA = '')
gpd_csv <- tail(list.files(path = "outputs/", 
           pattern = ".{10}_gpdpkglist-final",
           full.names = TRUE), 1)
if(length(gpd_csv) > 0){
  gpd_checks <- read.csv(file = gpd_csv)
} else{
  stop("csv file not found")
}
gpd_checks |> 
  dplyr::filter(!is.na(fun)) |>
  dplyr::select(package,
                  location,
                  density.comment, 
                  distribution.comment) |>
  dplyr::mutate(location = ifelse(is.na(location),
                                    "no",
                                    "yes")) |>
  dplyr::arrange(tolower(package)) |>
  knitr::kable(col.names = c("package",
                                  "location",
                                  "density",
                                  "distribution function"),
               row.names = FALSE,
               booktabs = TRUE) |>
  kableExtra::kable_styling()


#| label: tbl-gev1
#| tbl-cap: "Evaluation of  generalized extreme value density and distribution functions."
#| eval: true
#| echo: false
#| message: false
#| warning: false
#| cache: true
#| error: false

options(knitr.kable.NA = '')
gev_csv <- tail(list.files(path = "outputs/", 
           pattern = ".{10}_gevpkglist-final.csv", full.names = TRUE), 1)
gev_checks <- read.csv(file = gev_csv)
gev_checks |> 
  dplyr::filter(!is.na(fun)) |>
  dplyr::select(package,
                  density.comment, 
                  distribution.comment) |>
   dplyr::arrange(tolower(package)) |>
  knitr::kable(col.names = c("package",
                             "density",
                             "distribution function"),
               row.names = FALSE,
               booktabs = TRUE) |>
  kableExtra::kable_styling()


#| label: fig-gpfit-comput
#| cache: true
#| echo: false
#| eval: true
#| message: false
#| warning: false
# library(dplyr::, warn.conflicts = FALSE)
for(file in list.files(path = "outputs/", 
                       pattern = ".{11}GPD-fit.{3}.RData", 
                       full.names = TRUE)){
   load(file)
}

ids <- with(res_neg,
            which(apply(fullres[, "n20", "shape", ], 2, function(x) {
              sum(is.na(x))
            }) == dim(fullres)[1]))

# With 400 observations and a threshold set to the
# theoretical 95% percentile, we have on average
# 20 observations for the fit.
results <- res_neg$fullres[, "n50", , -ids]

g1 <- results[,  "gradshape", ] |>
  tibble::as_tibble() |>
  tidyr::pivot_longer(cols = tidyselect::everything(),
               names_to = "package",
               values_to = "gradient") |>
  dplyr::mutate(boundary = 
                  factor(rep(results[,  "shape", "mev"] < -0.9999, each = length(unique(package)))),
                gradient = ifelse(is.infinite(gradient), NA, gradient)) |>
  dplyr::arrange(tolower(package)) |>
  ggplot(mapping = aes(y = forcats::fct_rev(package),
                       x = abs(gradient),
                       fill = boundary)) +
  scale_fill_grey() +
  ggdist::stat_slab(normalize = "xy",
                    n = 1001,) +
  theme_minimal() +
  theme(legend.position = "none") +
   labs(x = "absolute gradient, shape (log scale)",
        y = "") +
   scale_x_log10(labels = scales::label_log(digits = 2),
                 limits = c(1e-8, 1e2),
                 oob = scales::squish,
                 breaks = c(1e-8, 1e-6, 1e-7, 1e-5, 1e-4, 1e-3,  1e-2, 1e-1, 1e0, 1e1),
                 minor_breaks = c(1e-10,1e-9,1e-8,1e-7,1e-6, 1e-5, 1e-4, 1e-3, 1e-2, 1e-1, 1e0, 1e1, 1e2)) +
   geom_vline(xintercept = 1e-3)

# Find which is the maximum likelihood
shape_mle <- llmax <- scale_mle <- vector("numeric", 1000L)
wmax <- vector("integer", 1000L)
for (i in seq_along(shape_mle)) {
  wmax[i] <- which.max(results[i, "loglik", ])
  shape_mle[i] <- results[i, "shape", wmax[i]]
  scale_mle[i] <- results[i, "scale", wmax[i]]
  llmax[i] <- results[i, "loglik", wmax[i]]
}


g2 <- tibble::as_tibble(results[, "loglik", ] - llmax) |>
  dplyr::select(names(which(apply(abs(results[, "loglik", ] - llmax), 2, quantile, 0.9, na.rm = TRUE) > 1e-4))) |>
  tidyr::pivot_longer(cols = tidyselect::everything(),
               names_to = "package",
               values_to = "loglik") |>
  dplyr::arrange(tolower(package)) |>
  dplyr::filter(rep(results[, "shape", "mev"] >= -0.9999, each = length(unique(package)))) |>
  dplyr::mutate(loglik = ifelse(is.finite(loglik), loglik, NA)) |>
  ggplot(mapping = aes(x = abs(loglik),
                       y = forcats::fct_rev(package))) +
  ggdist::stat_slab(slab_type = "ccdf") +
  theme_classic() +
  labs(x = "difference to maximum log likelihood",
       y = "")


g3 <- tibble::as_tibble(results[,  "shape", ]) |>
  tidyr::pivot_longer(cols = tidyselect::everything(),
               names_to = "package",
               values_to = "shape") |>
  dplyr::arrange(tolower(package)) |>
  ggplot(mapping = aes(x = shape, y = forcats::fct_rev(package))) +
  ggdist::stat_dots() +
  theme_classic() +
  labs(x = "shape",
       y = "")
results <- res_exp$fullres[, "n50", , -ids]
# Find which is the maximum likelihood
shape_mle <- llmax <- scale_mle <- vector("numeric", 1000L)
wmax <- vector("integer", 1000L)
for (i in seq_along(shape_mle)) {
  wmax[i] <- which.max(results[i, "loglik", ])
  shape_mle[i] <- results[i, "shape", wmax[i]]
  scale_mle[i] <- results[i, "scale", wmax[i]]
  llmax[i] <- results[i, "loglik", wmax[i]]
}

g4 <- tibble::as_tibble(results[, "loglik", ] - llmax) |>
  dplyr::select(names(which(apply(abs(results[, "loglik", ] - llmax), 2, quantile, 0.9, na.rm = TRUE) > 1e-4))) |>
  tidyr::pivot_longer(cols = tidyselect::everything(),
               names_to = "package",
               values_to = "loglik") |>
  dplyr::arrange(tolower(package)) |>
  dplyr::filter(rep(results[, "shape", "mev"] >= -0.9999, each = length(unique(package)))) |>
  dplyr::mutate(loglik = ifelse(is.finite(loglik), loglik, NA)) |>
  ggplot(mapping = aes(x = abs(loglik),
                       y = forcats::fct_rev(package))) +
  ggdist::stat_slab(slab_type = "ccdf") +
  theme_classic() +
  labs(x = "difference to maximum log likelihood",
       y = "")


g5 <- tibble::as_tibble(results[,  "shape", ]) |>
  tidyr::pivot_longer(cols = tidyselect::everything(),
               names_to = "package",
               values_to = "shape") |>
  dplyr::arrange(tolower(package)) |>
  ggplot(mapping = aes(x = shape, 
                       y = forcats::fct_rev(package))) +
  ggdist::stat_dots() +
  theme_classic() +
  labs(x = "shape",
       y = "")


#| label: fig-gpfit-grad
#| fig-cap: "Magnitude of the shape component of the score vector at the value returned by the optimization routine. The density plots are based on 1000 samples simulated from a generalized Pareto distribution with shape $\\xi =-0.5$ and scale $\\sigma=1000$, split by simulations yielding a boundary case ($\\widehat{\\xi} = -1$, gray) and regular case ($\\widehat{\\xi} > -1$, black); the $y$-axis scale for each package is different to ease visualization. Results for samples for which the numerical routines failed to converge or the gradient is unevaluated are not shown."
#| cache: false
#| echo: false
#| eval: true
#| message: false
#| warning: false
#| fig-align: 'center'
#| out-width: '80%'
#| fig-width: 8
#| fig-height: 6
g1


#| label: fig-gpfit-diffmle
#| fig-cap: "Difference between likelihood evaluated at parameters returned by routine and the maximum likelihood over all routines for generalized Pareto samples with negative shape (left) and exponential samples (right), both with large scale parameter $\\sigma=1000$. Results for samples for which the numerical routines failed to converge are not shown. Only packages with 90% percentile giving a discrepancy larger than $10^{-4}$ are shown."
#| cache: false
#| echo: false
#| eval: true
#| message: false
#| warning: false
#| fig-align: 'center'
#| out-width: '90%'
#| fig-width: 8
#| fig-height: 4
g2 + g4


#| label: fig-gpfit-shapepar
#| fig-cap: "Dot plots of shape parameters returned by optimization routine for generalized Pareto samples with negative shape (left) and exponential samples (right). Results for samples for which the numerical routines failed to converge are not shown."
#| cache: false
#| echo: false
#| eval: true
#| message: false
#| warning: false
#| fig-align: 'center'
#| out-width: '90%'
#| fig-width: 8
#| fig-height: 4
g3 + g5


#| label: tab-nonconv
#| tab-cap: "Number of failures (out of 1000 simulations) per package for bounded ($\\xi=-0.5$), light-tailed ($\\xi=0$) and heavy-tailed ($\\xi=0.5$) samples from a generalized Pareto distribution."
#| cache: false
#| echo: false
#| eval: true
#| message: false
#| warning: false
nonconv_neg <- res_neg$noncvg[,1,]
nonconv_neg <- nonconv_neg[,colSums(nonconv_neg)>0] |>
  tibble::as_tibble(rownames = "n") |>
  dplyr::mutate(n = factor(n, 
                    levels = c("n20","n50","n100","n1000"),
                     labels = c(20L,50L, 100L, 1000L))) |>
  tidyr::pivot_longer(cols = -1,
                      names_to = "package",
                      values_to = "nfail")

nonconv_exp <- res_exp$noncvg[,1,]
nonconv_exp <- nonconv_exp[,colSums(nonconv_exp)>0] |>
  tibble::as_tibble(rownames = "n") |>
  dplyr::mutate(n = factor(n, 
                    levels = c("n20","n50","n100","n1000"),
                     labels = c(20L,50L, 100L, 1000L))) |>
  tidyr::pivot_longer(cols = -1,
                      names_to = "package",
                      values_to = "nfail")

nonconv_pos <- res_pos$noncvg[,1,]
nonconv_pos <- nonconv_pos[,colSums(nonconv_pos)>0] |>
  tibble::as_tibble(rownames = "n") |>
  dplyr::mutate(n = factor(n, 
                    levels = c("n20","n50","n100","n1000"),
                     labels = c(20L,50L, 100L, 1000L))) |>
  tidyr::pivot_longer(cols = -1,
                      names_to = "package",
                      values_to = "nfail")
fail_df <- dplyr::bind_rows(.id = "shape", 
                            neg = nonconv_neg, 
                            zero = nonconv_exp, 
                            pos = nonconv_pos)
# knitr::kable(fail_df, row.names = "n", col.names = c("shape","package"))
fail_qrm <- as.integer(unlist((fail_df |> dplyr::filter(shape == "neg", package == "QRM") |> dplyr::select(nfail))))
fail_evir <- as.integer(fail_df |> dplyr::filter(shape == "neg", n == "20", package == "evir") |> dplyr::select(nfail)) 


#| label: gevfit-computation
#| cache: true
#| echo: false
#| eval: true
#| message: false
#| warning: false
for(file in list.files(path = "outputs/", 
                       pattern = ".{11}GEV-fit-.{3}.RData",
                       full.names = TRUE)){
   load(file)
}
results <- res_neg[, "n20",, ]

shape_mle <- llmax <- scale_mle <- vector("numeric", 1000L)
wmax <- vector("integer", 1000L)
for (i in seq_along(shape_mle)) {
  wmax[i] <- which.max(results[i, "loglik", results[i, "shape", ] > -1])
  shape_mle[i] <- results[i, "shape", wmax[i]]
  scale_mle[i] <- results[i, "scale", wmax[i]]
  llmax[i] <- results[i, "loglik", wmax[i]]
}


g6 <- tibble::as_tibble(results[, "loglik", ] - llmax) |>
  dplyr::select(names(which(apply(abs(results[, "loglik", ] - llmax), 2, quantile, 0.9, na.rm = TRUE) > 1e-4))) |>
  tidyr::pivot_longer(cols = tidyselect::everything(),
               names_to = "package",
               values_to = "loglik") |>
  dplyr::arrange(tolower(package)) |>
  dplyr::filter(rep(results[, "shape", "mev"] >= -0.9999, each = length(unique(package)))) |>
  dplyr::mutate(loglik = ifelse(is.finite(loglik), loglik, NA)) |>
  ggplot(mapping = aes(x = abs(loglik),
                       y = forcats::fct_rev(package))) +
  ggdist::stat_slab(slab_type = "ccdf", limits = c(0, 8)) +
  theme_classic() +
  labs(x = "difference to maximum log likelihood",
       y = "") +
  scale_x_continuous(limits = c(0, 8))
g7 <- tibble::as_tibble(results[,  "shape", ]) |>
  tidyr::pivot_longer(cols = tidyselect::everything(),
               names_to = "package",
               values_to = "shape") |>
  dplyr::arrange(tolower(package)) |>
  ggplot(mapping = aes(x = shape, y = package)) +
  ggdist::stat_dots() +
  theme_classic() +
  labs(x = "shape",
       y = "") +
  scale_x_continuous(limits = c(-2,1),oob = scales::squish)


results <- res_zer[, "n20", ,]
# Find which is the maximum likelihood
shape_mle <- llmax <- scale_mle <- vector("numeric", 1000L)
wmax <- vector("integer", 1000L)
for (i in seq_along(shape_mle)) {
  wmax[i] <- which.max(results[i, "loglik", ])
  shape_mle[i] <- results[i, "shape", wmax[i]]
  scale_mle[i] <- results[i, "scale", wmax[i]]
  llmax[i] <- results[i, "loglik", wmax[i]]
}

g8 <- tibble::as_tibble(results[, "loglik", ] - llmax) |>
  dplyr::select(names(which(apply(abs(results[, "loglik", ] - llmax), 2, quantile, 0.9, na.rm = TRUE) > 1e-4))) |>
  tidyr::pivot_longer(cols = tidyselect::everything(),
               names_to = "package",
               values_to = "loglik") |>
  dplyr::arrange(tolower(package)) |>
  dplyr::filter(rep(results[, "shape", "mev"] >= -0.9999, each = length(unique(package)))) |>
  dplyr::mutate(loglik = ifelse(is.finite(loglik), loglik, NA)) |>
  ggplot(mapping = aes(x = abs(loglik),
                       y = forcats::fct_rev(package))) +
  ggdist::stat_slab(slab_type = "ccdf") +
  theme_classic() +
  labs(x = "difference to maximum log likelihood",
       y = "")


g9 <- tibble::as_tibble(results[,  "shape", ]) |>
  tidyr::pivot_longer(cols = tidyselect::everything(),
               names_to = "package",
               values_to = "shape") |>
  dplyr::arrange(tolower(package)) |>
  # filter(package != "tea") |> # outlier
  ggplot(mapping = aes(x = shape, y = package)) +
  ggdist::stat_dots() +
  # geom_boxplot() +
  # scale_x_continuous(limits = c(-1, 2)) +
  # ggdist::stat_halfeye() +
  theme_classic() +
  labs(x = "shape",
       y = "") +
  scale_x_continuous(limits = c(-2,1),
                     oob = scales::squish)


#| label: fig-gevfit-diffmle
#| fig-cap: "Difference between likelihood evaluated at parameters returned by routine and the maximum likelihood over all routines for generalized extreme value samples with negative shape (left) and exponential samples (right), both with large scale parameter $\\sigma=1000$. Results for samples for which the numerical routines failed to converge are not shown.  Only packages with 90% percentile giving a discrepancy larger than $10^{-4}$ are shown."
#| cache: false
#| echo: false
#| eval: true
#| message: false
#| warning: false
#| fig-align: 'center'
#| out-width: '90%'
#| fig-width: 8
#| fig-height: 4
g6 + g8


#| label: fig-gevfit-shapepar
#| fig-cap: "Dot plots of shape parameters returned by optimization routine for generalized extreme value samples with negative shape (left) and Gumbel samples (right). Results for samples for which the numerical routines failed to converge are not shown."
#| cache: false
#| echo: false
#| eval: true
#| message: false
#| warning: false
#| fig-align: 'center'
#| out-width: '90%'
#| fig-width: 8
#| fig-height: 4
g7 + g9


#| label: tbl-gevnonconv
#| tbl-cap: "Number of failures for the optimization routine for maximum likelihood-based estimation of the generalized extreme value model (out of 1000 simulations)."
#| tbl-subcap:
#|   - "bounded tail"
#|   - "light tail"
#|   - "heavy tail"
#| layout-ncol: 1
#| cache: true
#| eval: true
#| echo: false
#| message: false
#| error: false
for(file in list.files(path = "outputs/", 
                       pattern = ".{11}GEV-fit-.{3}.RData",
                       full.names = TRUE)){
   load(file)
}
nonconv_neg <- apply(res_neg[,, "shape",], 2:3, function(x){sum(is.na(x))})
nonconv_neg <- nonconv_neg[,colSums(nonconv_neg)>0] |>
  tibble::as_tibble(rownames = "n") |>
  dplyr::mutate(n = factor(n, 
                    levels = c("n20","n50","n100","n1000"),
                     labels = c(20L,50L, 100L, 1000L))) #|>
  # tidyr::pivot_longer(cols = -1,
  #                     names_to = "package",
  #                     values_to = "nfail")
knitr::kable(nonconv_neg,
             booktabs = TRUE) |>
  kableExtra::kable_styling()
nonconv_zer <- apply(res_zer[,, "shape",], 2:3, function(x){sum(is.na(x))})
nonconv_zer <- nonconv_zer[,colSums(nonconv_zer)>0] |>
  tibble::as_tibble(rownames = "n") |>
  dplyr::mutate(n = factor(n, 
                    levels = c("n20","n50","n100","n1000"),
                     labels = c(20L,50L, 100L, 1000L))) #|>
  # tidyr::pivot_longer(cols = -1,
  #                     names_to = "package",
  #                     values_to = "nfail")
knitr::kable(nonconv_zer,
             booktabs = TRUE) |>
  kableExtra::kable_styling()
nonconv_pos <- apply(res_pos[,, "shape",], 2:3, function(x){sum(is.na(x))})
nonconv_pos <- nonconv_pos[,colSums(nonconv_pos)>0] |>
  tibble::as_tibble(rownames = "n") |>
  dplyr::mutate(n = factor(n, 
                    levels = c("n20","n50","n100","n1000"),
                     labels = c(20L,50L, 100L, 1000L))) #|>
#   tidyr::pivot_longer(cols = -1,
#                       names_to = "package",
#                       values_to = "nfail")
# fail_df <- dplyr::bind_rows(.id = "shape", 
#                             neg = nonconv_neg, 
#                             zero = nonconv_zer, 
#                             pos = nonconv_pos)
knitr::kable(nonconv_pos,
             booktabs = TRUE) |>
  kableExtra::kable_styling()

