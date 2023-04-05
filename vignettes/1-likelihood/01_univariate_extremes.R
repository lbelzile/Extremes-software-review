#| label: setup
#| echo: false
#| eval: true
library(ggplot2)
library(patchwork)


#| label: tbl-gp1
#| tbl-cap: "Evaluation of density and distribution function for the generalized Pareto model."
#| eval: true
#| echo: false
#| message: false
#| warning: false
#| cache: true
setwd(here::here())
options(knitr.kable.NA = '')
gpd_csv <- "outputs/2022_03_01_gpdpkglist-final.csv"
gpd_checks <- read.csv(file = gpd_csv)
gpd_checks |> 
  dplyrfilter(!is.na(fun)) |>
  dplyrselect(package,
                  location,
                  density.comment, 
                  distribution.comment) |>
  dplyrmutate(location = ifelse(is.na(location),
                                    "no",
                                    "yes")) |>
   dplyr::arrange("package") |>
  knitr::kable(col.names = c("package",
                                  "location",
                                  "density",
                                  "distribution function"),
                    row.names = FALSE)


#| label: tbl-gev1
#| tbl-cap: "Evaluation of density and distribution function for the generalized extreme value distribuion."
#| eval: true
#| echo: false
#| message: false
#| warning: false
#| cache: true
#| error: false

options(knitr.kable.NA = '')
gev_csv <- "outputs/2022_03_01_gevpkglist-final.csv"
gev_checks <- read.csv(file = gev_csv)
gev_checks |> 
  dplyrfilter(!is.na(fun)) |>
  dplyrselect(package,
                  density.comment, 
                  distribution.comment) |>
   dplyr::arrange("package") |>
  knitr::kable(col.names = c("package",
                             "density",
                             "distribution function"),
                    row.names = FALSE)


#| label: fig-gpfit-comput
#| cache: true
#| echo: false
#| eval: true
#| message: false
#| warning: false
# library(dplyr, warn.conflicts = FALSE)
for(file in list.files(path = "outputs/", 
                       pattern = ".{11}GPD-fit.{3}.RData")){
   load(paste0("outputs/",file))
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
                ifelse(is.infinite(gradient), NA, gradient)) |>
  # filter(boundary == "no") |>
   dplyr::arrange("package") |>
  ggplot(mapping = aes(y = package,
                       x = log(abs(gradient)),
                       fill = boundary)) +
  scale_fill_manual(values = MetBrewer::met.brewer("OKeeffe1",
                                                   type = "discrete",
                                                   n = 2)) +
  ggdist::stat_slab(normalize = "xy", 
                    limits  = c(-30,30)) +
  theme_minimal() +
  theme(legend.position = "none") +
  labs(x = "log absolute gradient (shape)",
       y = "") +
  geom_vline(xintercept = log(0.001))

g2 <- results[,  "gradscale", ] |>
  tibble::as_tibble() |>
  tidyr::pivot_longer(cols = tidyselect::everything(),
               names_to = "package",
               values_to = "gradient") |>
  dplyr::mutate(boundary = factor(rep(results[, "shape", "mev"] < -0.9999, each = length(unique(package)))),
                package = factor(package)) |>
  # filter(boundary == "no") |>
  ggplot(mapping = aes(y = package,
                       x = log(abs(gradient)),
                       fill = boundary)) +
  scale_fill_manual(values = MetBrewer::met.brewer("OKeeffe1",
                                                   type = "discrete",
                                                   n = 2)) +
  ggdist::stat_slab(normalize = "xy", limits = c(-30,30)) +
  theme_minimal() +
  theme(legend.position = "none") +
  labs(x = "log absolute gradient (scale)",
       y = "") +
  geom_vline(xintercept = log(0.001))


# Find which is the maximum likelihood
shape_mle <- llmax <- scale_mle <- vector("numeric", 1000L)
wmax <- vector("integer", 1000L)
for (i in seq_along(shape_mle)) {
  wmax[i] <- which.max(results[i, "loglik", ])
  shape_mle[i] <- results[i, "shape", wmax[i]]
  scale_mle[i] <- results[i, "scale", wmax[i]]
  llmax[i] <- results[i, "loglik", wmax[i]]
}


g3a <- tibble::as_tibble(results[, "loglik", ] - llmax) |>
  tidyr::pivot_longer(cols = tidyselect::everything(),
               names_to = "package",
               values_to = "loglik") |>
  dplyr::filter(rep(results[, "shape", "mev"] >= -0.9999, each = length(unique(package)))) |>
  # filter(package != "tea") |> # outlier
  ggplot(mapping = aes(x = abs(loglik), y = package)) +
  geom_boxplot() +
  # ggdist::stat_slab() +
  theme_classic() +
  labs(x = "log likelihood difference to MLE",
       y = "")


g4a <- tibble::as_tibble(results[,  "shape", ]) |>
  tidyr::pivot_longer(cols = tidyselect::everything(),
               names_to = "package",
               values_to = "shape") |>
  # filter(package != "tea") |> # outlier
  ggplot(mapping = aes(x = shape, y = package)) +
  ggdist::stat_dots() +
  # geom_boxplot() +
  # scale_x_continuous(limits = c(-1, 2)) +
  # ggdist::stat_halfeye() +
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

g3b <- tibble::as_tibble(results[, "loglik", ] - llmax) |>
  tidyr::pivot_longer(cols = tidyselect::everything(),
               names_to = "package",
               values_to = "loglik") |>
  dplyr::filter(rep(results[, "shape", "mev"] >= -0.9999, each = length(unique(package)))) |>
  # filter(package != "tea") |> # outlier
  ggplot(mapping = aes(x = abs(loglik), y = package)) +
  geom_boxplot() +
  # ggdist::stat_slab() +
  theme_classic() +
  labs(x = "log likelihood difference to MLE",
       y = "")


g4b <- tibble::as_tibble(results[,  "shape", ]) |>
  tidyr::pivot_longer(cols = tidyselect::everything(),
               names_to = "package",
               values_to = "shape") |>
  # filter(package != "tea") |> # outlier
  ggplot(mapping = aes(x = shape, y = package)) +
  ggdist::stat_dots() +
  # geom_boxplot() +
  # scale_x_continuous(limits = c(-1, 2)) +
  # ggdist::stat_halfeye() +
  theme_classic() +
  labs(x = "shape",
       y = "")


#| label: fig-gpfit-grad
#| fig-cap: "Magnitude of the score vector at the value returned by the optimization routine. The density plots are based on 1000 samples simulated from a generalized Pareto distribution with shape $\\xi =-0.5$ and scale $\\sigma=1000$, split by simulations yielding a boundary case ($\\widehat{\\xi} = -1$, blue) and regular case ($\\widehat{\\xi} > -1$, red); the $y$-axis scale for each package is different to ease visualization. Results for samples for which the numerical routines failed to converge or the gradient is unevaluated are not shown. The vertical line indicates a tolerance of $0.001$"
#| cache: false
#| echo: false
#| eval: true
#| message: false
#| warning: false
g1 + g2


#| label: fig-gpfit-diffmle
#| fig-cap: "Difference between likelihood evaluated at parameters returned by routine and the maximum likelihood over all routines for generalized Pareto samples with negative shape (left) and exponential samples (right), both with large scale parameter $\\sigma=1000$. Results for samples for which the numerical routines failed to converge are not shown."
#| cache: false
#| echo: false
#| eval: true
#| message: false
#| warning: false
g3a + g3b


#| label: fig-gpfit-shapepar
#| fig-cap: "Dot plots of shape parameters returned by optimization routine for generalized Pareto samples with negative shape (left) and exponential samples (right). Results for samples for which the numerical routines failed to converge are not shown."
#| cache: false
#| echo: false
#| eval: true
#| message: false
#| warning: false
g4a + g4b


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
                       pattern = ".{11}GEV-fit-.{3}.RData")){
   load(paste0("outputs/",file))
}
results <- res_neg[, "n20",, ]
# 
# g1 <- results[,  "gradshape", ] |>
#   tibble::as_tibble() |>
#   tidyr::pivot_longer(cols = tidyselect::everything(),
#                names_to = "package",
#                values_to = "gradient") |>
#   dplyr::mutate(boundary = 
#                   factor(rep(results[,  "shape", "mev"] < -0.9999, each = length(unique(package)))),
#                 gradient = ifelse(is.infinite(gradient) | abs(gradient) > 1e10, NA, gradient)) |>
#   # filter(boundary == "no") |>
#   ggplot(mapping = aes(y = package,
#                        x = log(abs(gradient)),
#                        fill = boundary)) +
#   scale_fill_manual(values = MetBrewer::met.brewer("OKeeffe1",
#                                                    type = "discrete",
#                                                    n = 2)) +
#   ggdist::stat_slab(normalize = "xy", 
#                     limits  = c(-30,30)) +
#   theme_minimal() +
#   theme(legend.position = "none") +
#   labs(x = "log absolute gradient (shape)",
#        y = "") +
#   geom_vline(xintercept = log(0.001))+
#   scale_x_continuous(limits = c(-10,10))
# 
# g2 <- results[,  "gradscale", ] |>
#   tibble::as_tibble() |>
#   tidyr::pivot_longer(cols = tidyselect::everything(),
#                names_to = "package",
#                values_to = "gradient") |>
#   dplyr::mutate(boundary = factor(rep(results[, "shape", "mev"] < -0.9999, each = length(unique(package)))),
#                 package = factor(package)) |>
#   # filter(boundary == "no") |>
#   ggplot(mapping = aes(y = package,
#                        x = log(abs(gradient)),
#                        fill = boundary)) +
#   scale_fill_manual(values = MetBrewer::met.brewer("OKeeffe1",
#                                                    type = "discrete",
#                                                    n = 2)) +
#   ggdist::stat_slab(normalize = "xy", limits = c(-30,30)) +
#   theme_minimal() +
#   theme(legend.position = "none") +
#   labs(x = "log absolute gradient (scale)",
#        y = "") +
#   geom_vline(xintercept = log(0.001)) +
#   scale_y_continous(limits = c(-10,5))


# Find which is the maximum likelihood
shape_mle <- llmax <- scale_mle <- vector("numeric", 1000L)
wmax <- vector("integer", 1000L)
for (i in seq_along(shape_mle)) {
  wmax[i] <- which.max(results[i, "loglik", results[i, "shape", ] > -1])
  shape_mle[i] <- results[i, "shape", wmax[i]]
  scale_mle[i] <- results[i, "scale", wmax[i]]
  llmax[i] <- results[i, "loglik", wmax[i]]
}


g3a <- tibble::as_tibble(results[, "loglik", ] - llmax) |>
  tidyr::pivot_longer(cols = tidyselect::everything(),
               names_to = "package",
               values_to = "loglik") |>
  dplyr::filter(rep(results[, "shape", "mev"] >= -0.9999, each = length(unique(package)))) |>
  ggplot(mapping = aes(x = abs(loglik), y = package)) +
  geom_boxplot() +
  # ggdist::stat_slab() +
  theme_classic() +
  labs(x = "log likelihood difference to MLE",
       y = "") + 
  scale_x_continuous(limits = c(0, 15))
g4a <- tibble::as_tibble(results[,  "shape", ]) |>
  tidyr::pivot_longer(cols = tidyselect::everything(),
               names_to = "package",
               values_to = "shape") |>
  # filter(package != "tea") |> # outlier
  ggplot(mapping = aes(x = shape, y = package)) +
  ggdist::stat_dots() +
  # geom_boxplot() +
  # scale_x_continuous(limits = c(-1, 2)) +
  # ggdist::stat_halfeye() +
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

g3b <- tibble::as_tibble(results[, "loglik", ] - llmax) |>
  tidyr::pivot_longer(cols = tidyselect::everything(),
               names_to = "package",
               values_to = "loglik") |>
  dplyr::filter(rep(results[, "shape", "mev"] >= -0.9999, each = length(unique(package)))) |>
  # filter(package != "tea") |> # outlier
  ggplot(mapping = aes(x = abs(loglik), y = package)) +
  geom_boxplot() +
  # ggdist::stat_slab() +
  theme_classic() +
  labs(x = "log likelihood difference to MLE",
       y = "") +
  scale_x_continuous(limits = c(0,15))


g4b <- tibble::as_tibble(results[,  "shape", ]) |>
  tidyr::pivot_longer(cols = tidyselect::everything(),
               names_to = "package",
               values_to = "shape") |>
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
#| fig-cap: "Difference between likelihood evaluated at parameters returned by routine and the maximum likelihood over all routines for generalized extreme value samples with negative shape (left) and exponential samples (right), both with large scale parameter $\\sigma=1000$. Results for samples for which the numerical routines failed to converge are not shown."
#| cache: false
#| echo: false
#| eval: true
#| message: false
#| warning: false
g3a + g3b


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
                       pattern = ".{11}GEV-fit-.{3}.RData")){
   load(paste0("outputs/",file))
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
knitr::kable(nonconv_neg)
nonconv_zer <- apply(res_zer[,, "shape",], 2:3, function(x){sum(is.na(x))})
nonconv_zer <- nonconv_zer[,colSums(nonconv_zer)>0] |>
  tibble::as_tibble(rownames = "n") |>
  dplyr::mutate(n = factor(n, 
                    levels = c("n20","n50","n100","n1000"),
                     labels = c(20L,50L, 100L, 1000L))) #|>
  # tidyr::pivot_longer(cols = -1,
  #                     names_to = "package",
  #                     values_to = "nfail")
knitr::kable(nonconv_zer)
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
knitr::kable(nonconv_pos)

