library(ggplot2)
library(patchwork)
library(tidyverse)

setwd(here::here())
# Load results from simulation study
# Fitting generalized Pareto to 1000 
# simulated Gamma samples of size 200 to 9000
load("Simulation_gamma.RData")

ids <- with(res1gamma,
            which(apply(fullres[, "n500", "shape", ], 2, function(x) {
               sum(is.na(x))
            }) == dim(fullres)[1]))
results <- res1gamma$fullres[, "n500", , -ids]


g1 <- tibble::as_tibble(results[,  "shape", ]) |>
   tidyr::pivot_longer(cols = tidyselect::everything(),
                       names_to = "package",
                       values_to = "shape") |>
   # filter(package != "tea") |> # outlier
   ggplot(mapping = aes(x = shape, y = forcats::fct_rev(package))) +
   ggdist::stat_dots(normalize = "xy") +
   # geom_boxplot() +
   # scale_x_continuous(limits = c(-1, 2)) +
   # ggdist::stat_halfeye() +
   theme_classic() +
   labs(x = "shape",
        y = "")
        
g2 <- results[,  "gradshape", ] |>
   tibble::as_tibble() |>
   tidyr::pivot_longer(cols = tidyselect::everything(),
                       names_to = "package",
                       values_to = "gradient") |>
   dplyr::mutate(boundary =
                    factor(rep(results[,  "shape", "mev"] < -0.9999, each = length(unique(package)))),
                 gradient = ifelse(is.infinite(gradient), NA, gradient)) |>
   dplyr::filter(boundary == FALSE) |>
   ggplot(mapping = aes(y = forcats::fct_rev(package),
                        # fill = boundary,
                        x = abs(gradient))
          ) +
   scale_fill_manual(values = c("grey10","grey70")) +
   ggdist::stat_slab(n = 2001) +
   theme_minimal() +
   theme(legend.position = "none") +
   labs(x = "absolute gradient, shape (log scale)",
        y = "") +
   scale_x_log10(labels = scales::label_log(digits = 2),
                 limits = c(1e-8, 1e2),
                 oob = scales::squish,
                 breaks = c(1e-8, 1e-6, 1e-7, 1e-5, 1e-4, 1e-3,  1e-2, 1e-1, 1e0, 1e1),
                 minor_breaks = c(1e-10,1e-9,1e-8,1e-7,1e-6, 1e-5, 1e-4, 1e-3, 1e-2, 1e-1, 1e0, 1e1, 1e2)) +
   geom_vline(xintercept = 1e-3) +
  theme_classic()

# Find which is the maximum likelihood
shape_mle <- llmax <- scale_mle <- vector("numeric", 1000L)
wmax <- vector("integer", 1000L)
for (i in seq_along(shape_mle)) {
   wmax[i] <- which.max(results[i, "loglik", ])
   shape_mle[i] <- results[i, "shape", wmax[i]]
   scale_mle[i] <- results[i, "scale", wmax[i]]
   llmax[i] <- results[i, "loglik", wmax[i]]
}



t3 <- tibble::as_tibble(results[, "loglik", ] - llmax) |>
   tidyr::pivot_longer(cols = tidyselect::everything(),
                       names_to = "package",
                       values_to = "loglik") |>
   dplyr::filter(rep(results[, "shape", "mev"] >= -0.9999, each = length(unique(package)))) |>
   dplyr::group_by(package) |>
   dplyr::summarize(percentage = mean(abs(loglik) > 0.05))

pdf("fig2_GP-gamma_20exc.pdf",
    width = 10, height = 8)
g1 + g2
dev.off()

