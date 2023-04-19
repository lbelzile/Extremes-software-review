#| label: extremes_benchmark_setup
#| echo: false
#| message: false
#| error: false
#| warning: false
if (!"this.path" %in% installed.packages()[, "Package"]){
  install.packages("this.path")
}
setwd(this.path::here())
if(!dir.exists("outputs")){
 mkdir("outputs")
}
data(venice, package = "mev")
data(eskrain, package = "mev")
suppressPackageStartupMessages(library(cmdstanr))
# Compile C++ model
sm1 <-  cmdstan_model("gev.stan")
sm2 <-  cmdstan_model("gpd.stan")


#| label: extremes_benchmark1
#| echo: false
#| message: true
#| warning: false
#| error: false
set.seed(12345)
if(file.exists("outputs/benchmark1.RData")){
  load("outputs/benchmark1.RData")
} else{
  benchmark1 <- 
    bench::mark(iterations = 100,
              filter_gc = TRUE,
              check = FALSE,
  texmex = texmex::evm(
    y = r1, 
    mu = ~ scale(year),
    data = venice, 
    burn = 1000,
    iter  = 4.1e4,
    family = texmex::gev, 
    method = "sim",
    verbose = FALSE),
  extRemes = {
    prelim_opt <- extRemes::fevd(
      x = r1,
      data = venice[,1:2],
      location.fun = ~scale(year),
      type = "GEV", 
      use.phi = TRUE,
      method = "MLE")
  
  # Initial trials with the default parameter revealed
  # problems: while the model starts at the MLE (so close to the stationary), the default standard deviation of the normal random walk proposals are particulary ill-suited to this example.
  extRemes <- extRemes::fevd(
      x = r1,
      data = venice[,1:2],
      location.fun = ~scale(year),
      type = "GEV", 
      use.phi = TRUE,
      method = "Bayesian",
      iter = 4e4L+1e3L,
      priorParams = list(
        m = rep(0,4),
        v = c(1000, 1000, 1000, 1)
        ),
      proposalParams = list(
        mean = rep(0,4),
        sd = sqrt(diag(solve(prelim_opt$results$hessian)))
        )
  )
  },
  evdbayes = {
    prelim_opt <- with(venice, evd::fgev(
      x = r1, 
      nsloc = data.frame(x = scale(year)))
    )
    evdbayes::posterior(
    n = 4e4L+999,
    init = prelim_opt$estimate[c(1,3,4,2)],
    psd = prelim_opt$std.err[c(1,3,4,2)],
    lh = "gev",
    data = venice$r1,
    burn = 1000L,
    thin = 4,
    trend = scale(venice$year),
    prior = evdbayes::prior.norm(trendsd = 5,
                       mean = rep(0,3),
                       cov = diag(c(1000,1000,1)))
      
  )
  },
  stan ={ 
    prelim_opt <- with(venice, evd::fgev(
      x = r1, 
      nsloc = data.frame(x = scale(year)))
    )
    invisible(capture.output( 
    suppressWarnings(
      sm1$sample(data = list(N = nrow(venice),
                    p_loc = 2L,
                    priorvar = c(100,100,1),
                    p_scale = 1L,
                    y = venice$r1,
                    X_loc = cbind(rep(1, nrow(venice)),
                                  scale(venice$year)),
                    X_scale = cbind(rep(1, nrow(venice)))),
                  init = list(list(beta_loc = c(1.42011,1.35351),
                                   beta_scale = c(0.96582),
                                   xi = 0.04071)),
                   chains = 1,
                   parallel_chains = 1L,
                   iter_warmup = 500L,
                   iter_sampling = 1e4L,
                   refresh = 0L, # don't print progress
                   show_messages = FALSE))))
  },
  
  )
  save(benchmark1, file = "outputs/benchmark1.RData", version = 2)
}


#| label: extremes_benchmark2
#| echo: false
#| message: false
#| warnings: false
#| error: false
set.seed(12345)
if(file.exists("outputs/benchmark2.RData")){
  load("outputs/benchmark2.RData")
} else{
  benchmark2 <- 
    suppressWarnings(bench::mark(iterations = 100,
              filter_gc = FALSE,
              check = FALSE,
  texmex = texmex::evm(
    y = eskrain, 
    th = 30, 
    family = texmex::gpd, 
    method = "sim",
    verbose = FALSE),
  revdbayes = revdbayes::rpost_rcpp(
    n = 1e4L, 
    model = "gp", 
    data = eskrain, 
    thresh = 30, 
    prior = revdbayes::set_prior(prior = "norm", 
                                 model = "gp", 
                                 mean = rep(0,2),
                                 cov = diag(c(1000,1)))
  ),
  MCMC4Extremes = invisible(capture.output(
    mod2_MCMC4_s = MCMC4Extremes::gpdp(
      data = eskrain, 
      threshold = 30, 
      int = 1e4L))), # this prints the whole progression!
  
  extRemes = extRemes::fevd(
      x = as.vector(eskrain),
      threshold = 30,
      type = "GP", 
      use.phi = TRUE,
      method = "Bayesian",
      iter = 1e4L+500L),
  stan = {invisible(capture.output( 
    suppressWarnings(
    sm2$sample(data = list(
                             N = length(eskrain),
                             thresh = 30,
                             y = as.vector(eskrain)),
                   chains = 1,
                   parallel_chains = 1L,
                   iter_warmup = 500L,
                   iter_sampling = 1e4L,
                   refresh = 0L, # don't print progress
                   show_messages = FALSE)
  )))
  }
  ))
  save(benchmark2, file = "outputs/benchmark2.RData", version = 2)
}


#| label: fig-benchmark
#| fig.cap: "Swarm plot of speed (including preliminary optimization if necessary) of different numerical routines for a generalized extreme value model with linear trend in location fitted to the Venice sea level data (left) and the generalized Pareto distribution fitted to the Eskdalemuir rainfall data (right)."
#| eval: true
#| echo: false
#| error: false
#| message: false
#| warning: false
#| fig-width: 10
#| fig-height: 4
#| out-width: '90%'
#| fig-align: "center"
library(bench)
library(ggplot2)
library(patchwork)
 g1 <- autoplot(benchmark1, color = "black") +
   theme_classic() +
   labs(y = "time (log scale)",
        x = "package",
        subtitle = "Venice sea level data")
   
g2 <- autoplot(benchmark2, color = "black") + 
   theme_classic() +
   labs(y = "time (log scale)",
        x = "",
        subtitle = "Eskdalemuir rainfall data")
g1 + g2


#| label: evaless
#| echo: false
#| eval: true
#| message: false
#| warning: false
#| error: false
#| results: 'hide'
set.seed(12345)
if(file.exists("outputs/ESS.RData")){
  load("outputs/ESS.RData")
} else{
set.seed(12345)
mod1_texmex <- texmex::evm(
  y = r1, 
  mu = ~ scale(year),
  data = venice, 
  chains = 1,
  burn = 500,
  trace = Inf,
  verbose = FALSE,
  thin = 4,
  family = texmex::gev, 
  method = "sim")

# Estimate effective sample size
# Fit an AR process via glm
texmex_ess1 <- with(mod1_texmex,
coda::effectiveSize(coda::as.mcmc(param))/(thin*nrow(param))
)
## While we cannot know if we have converged to the target posterior distribution, the chains appear stationary
# par(mfrow = c(2,2))
# coda::traceplot(coda::as.mcmc(mod1_texmex$param))
# dev.off()

mod2_texmex = texmex::evm(
  y = eskrain, 
  th = 30, 
  chains = 1,
  family = texmex::gpd, 
  method = "sim",
  verbose = FALSE)
texmex_ess2 <- with(mod2_texmex,
coda::effectiveSize(coda::as.mcmc(param))/(thin*nrow(param))
)


invisible(capture.output(
  mod2_MCMC4 <- MCMC4Extremes::gpdp(
    data = eskrain, 
    threshold = 30, 
    int = 1e4L)))
MCMC4_ess2 <- with(mod2_MCMC4,
coda::effectiveSize(coda::as.mcmc(posterior))/(10*nrow(posterior)))

# Check traceplots for convergence
# par(mfrow = c(1,2))
# coda::traceplot(coda::as.mcmc(mod2_MCMC4$posterior))
# dev.off()
# 
# Fast considering the number of observations. But crude and inefficient.
# Includes a burnin of 50K observations! (int*thin / 2, contrary to what is in the documentation) - not customizable

prelim_opt <- extRemes::fevd(
    x = r1,
    data = venice[,1:2],
    location.fun = ~scale(year),
    type = "GEV", 
    use.phi = TRUE,
    method = "MLE")

# Initial trials with the default parameter revealed
# problems: while the model starts at the MLE (so close to the stationary), the default standard deviation of the normal random walk proposals are particulary ill-suited to this example.
mod1_extRemes <- extRemes::fevd(
    x = r1,
    data = venice[,1:2],
    location.fun = ~scale(year),
    type = "GEV", 
    use.phi = TRUE,
    method = "Bayesian",
    iter = 1e4L+500L,
    priorParams = list(
      m = rep(0,4),
      v = c(1000, 1000, 1000, 1)
      ),
    proposalParams = list(
      mean = rep(0,4),
      sd = sqrt(diag(solve(prelim_opt$results$hessian)))
      )
)
# Produce traceplots
# par(mfrow = c(2,4))
# plot(mod1_extRemes, type = "trace")
# dev.off()

extRemes_ess1 <- with(mod1_extRemes,
coda::effectiveSize(coda::as.mcmc(results[-(1:500),1:4]))/1e4)
# Effective sample size is very small...
# Trace plots reveal lack of stationarity
# with default tuning parameters
# 
# With adapted proposals (and vague priors), the output seems satisfactory, but the effective sample size is subpar with other methods.

mod2_extRemes <- extRemes::fevd(
    x = as.vector(eskrain),
    threshold = 30,
    type = "GP", 
    use.phi = TRUE,
    method = "Bayesian",
    iter = 1e4L+500L)
# Remove burnin
extRemes_ess2 <- with(mod2_extRemes,
coda::effectiveSize(
  coda::as.mcmc(results[-(1:500),1:2]))/1e4)

prelim_opt <- with(venice, evd::fgev(
    x = r1, 
    nsloc = data.frame(x = scale(year)))
  )
mod1_evdbayes <- evdbayes::posterior(
  n = 1e4L+500,
  init = prelim_opt$estimate[c(1,3,4,2)],
  psd = prelim_opt$std.err[c(1,3,4,2)],
  lh = "gev",
  data = venice$r1,
  burn = 499L,
  trend = scale(venice$year),
  prior = evdbayes::prior.norm(trendsd = 5,
                     mean = rep(0,3),
                     cov = diag(c(1000,1000,1)))
    
)
evdbayes_ess1 <- coda::effectiveSize(coda::as.mcmc(mod1_evdbayes))/nrow(mod1_evdbayes)
# Generalized Pareto model with evdbayes includes
# a location parameter, so not comparable...

mod1_stan <- sm1$sample(data = list(N = nrow(venice),
                  p_loc = 2L,
                  priorvar = c(100,100,1),
                  p_scale = 1L,
                  y = venice$r1,
                  X_loc = cbind(rep(1, nrow(venice)),
                                scale(venice$year)),
                  X_scale = cbind(rep(1, nrow(venice)))),
                init = list(list(beta_loc = c(1.42011,1.35351),
                                 beta_scale = c(0.96582),
                                 xi = 0.04071)),
                 chains = 1,
                 parallel_chains = 1L,
                 iter_warmup = 500L,
                 iter_sampling = 1e4L,
                 refresh = 0L, # don't print progress
                 show_messages = FALSE)

mod2_stan <- sm2$sample(data = list(
                           N = length(eskrain),
                           thresh = 30,
                           y = as.vector(eskrain)),
                 chains = 1,
                 parallel_chains = 1L,
                 iter_warmup = 500L,
                 iter_sampling = 10000L,
                 refresh = 0L, # don't print progress
                 show_messages = FALSE)
# Using the same method as before
stan_ess1 <- coda::effectiveSize(
  coda::as.mcmc(matrix(as.numeric(mod1_stan$draws(inc_warmup = FALSE)[,1,-1]), ncol = 4L)))/1e4

stan_ess2 <- coda::effectiveSize(
  coda::as.mcmc(matrix(as.numeric(mod2_stan$draws(inc_warmup = FALSE)[,1,-1]), ncol = 2L)))/1e4
# Since we can run independent chains, we could easily compute the effective sample size with multiple chains

save(stan_ess1, stan_ess2, evdbayes_ess1,
     extRemes_ess2, extRemes_ess1,
     MCMC4_ess2, texmex_ess1, texmex_ess2,
     file = "outputs/ESS.RData",
     version = 2)
}



#| label: tbl-ess
#| tbl-cap: "Effective sample size over number of iterations (percentage)."
#| tbl-subcap:
#|   - "nonstationary generalized extreme value model"
#|   - "generalized Pareto model"
#| layout-ncol: 1
#| echo: false
#| eval: true
ess1 <- data.frame(package = c("texmex",
                              "extRemes",
                              "evdbayes",
                              "STAN"))
ess1 <- cbind(ess1, 100*rbind(texmex_ess1,
                          extRemes_ess1,
                          evdbayes_ess1,
                          stan_ess1))
colnames(ess1)[-1] <- c("loc","loc (trend)","scale","shape")
ess2 <- data.frame(package = c("texmex",
                              "extRemes",
                              "STAN",
                              "MCMC4extremes"))
ess2 <- cbind(ess2, 100*rbind(texmex_ess2,
                          extRemes_ess2,
                          stan_ess2,
                          MCMC4_ess2))
colnames(ess2)[-1] <- c("scale","shape")

knitr::kable(ess1,
             digits = 1,
             align = c("lrrrr"),
             row.names = FALSE)

knitr::kable(ess2,
             digits = 1,
             align = c("lrr"),
             row.names = FALSE)


#| label: methods_bayesian
#| echo: false
#| warning: false
#| error: false
#| message: false
suppressPackageStartupMessages(
   library(texmex, 
           quietly = TRUE, 
           warn.conflicts = FALSE)
)
paste("Class 'evmSim' from the 'texmex' package")
methods(class = "evmSim")
detach("package:texmex", unload = TRUE)

 library(revdbayes, 
         quietly = TRUE, 
         warn.conflicts = FALSE)
paste("Class 'evpost' from the 'revdbayes' package")
methods(class = "evpost")
detach("package:revdbayes", unload = TRUE)

 suppressPackageStartupMessages(library(MCMC4Extremes, 
         quietly = TRUE, 
         warn.conflicts = FALSE))
paste("Classes 'gevp' and 'gpdp' from the 'MCMC4Extremes' package")
methods(class = "gevp")
methods(class = "gpdp")
 detach("package:MCMC4Extremes", unload = TRUE)
 library(extRemes, 
         quietly = TRUE, 
         warn.conflicts = FALSE)
paste("Class 'fecd'from the 'extRemes' package")
methods(class = "fevd")
paste("not all methods for this class are relevant.")
#not all methods for this class are pertinent for method = "Bayesian"
detach("package:extRemes", unload = TRUE)

