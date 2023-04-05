# Load libraries
library(ggplot2)
library(patchwork)
library(mev)
set.seed(1234)
dat <- mev::rgp(n = 10, shape = -0.8)
fitted <- mev::fit.gpd(dat, threshold = 0, show = TRUE)
# Empirical coefficient of variation
# Theoretical quantity defined as standard deviation/mean
sd(dat)/mean(dat)


prof_gpd_eta <- function(eta, xdat){
  -length(xdat) - sum(log(1-eta*xdat))-
    length(xdat)*log(-mean(log(1-eta*xdat))/eta)
}
grad_gpd_eta<- function(eta, xdat){
 sum(1/(1-eta*xdat))/eta - 
    length(xdat)/eta/(sum(log1p(-eta*xdat)))*(length(xdat) -
    sum(1/(1-eta*xdat)))
}

# Grid value for eta, excluding a neighborhood of zero
etagrid <- c(seq(-30, -1e-8, length = 201L),
             seq(1e-8, 1/max(dat)+1e-10,
                 length.out = 301L))
ll <- sapply(etagrid,
             FUN = prof_gpd_eta,
             xdat = dat)
# grad <- sapply(etagrid,
#                FUN = grad_gpd_eta,
#                xdat = dat)
# plot(etagrid, grad, ylim = c(-10,20), type = "l")
# For zero, we get xi=0 and exponential model,
# whose mle for scale is the sample mean
etagrid <- c(etagrid,0)

ll <- c(ll, sum(dexp(dat, rate = 1/mean(dat), log = TRUE)))

ll <- ll[order(etagrid)]
etagrid <- sort(etagrid)
xis <- sapply(etagrid, function(eta){mean(log(1-eta*dat))})
sub <- xis > -1
xis <- xis[sub]
etagrid <- etagrid[sub]
ll <- ll[sub]
# Plot the log likelihood
mllbound <- mev::gpd.ll(dat = dat, par = c(max(dat),-1))
mll <- pmax(max(ll, na.rm = TRUE), mllbound)

g1 <- ggplot(data = data.frame(eta = etagrid, cll = ll-mll),
       mapping = aes(x = eta, y = cll)) +
  geom_line()  +
  geom_point(data = data.frame(x = 1/max(dat), y = mllbound - mll),
             aes(x = x, y = y)) +
  theme_classic() +
  labs(
    x = expression(eta),
    y = "profile log likelihood")

dat <- evd::rgpd(n = 20, scale = 1, shape = 0.1)
etagrid <- c(seq(-4, -1e-8, length = 201L),
             seq(1e-8, 1/max(dat)+1e-10, length.out = 301L))
ll <- sapply(etagrid, FUN = prof_gpd_eta, xdat = dat)
# For zero, we get xi=0 and exponential model,
# whose mle for scale is the sample mean
etagrid <- c(etagrid,0)
ll <- c(ll, sum(dexp(dat, rate = 1/mean(dat), log = TRUE)))
ll <- ll[order(etagrid)]
etagrid <- sort(etagrid)
xis <- sapply(etagrid, function(eta){mean(log(1-eta*dat))})
sub <- xis > -1
xis <- xis[sub]
etagrid <- etagrid[sub]
ll <- ll[sub]
# Plot the log likelihood
mllbound <- mev::gpd.ll(dat = dat, par = c(max(dat),-1))
mll <- pmax(max(ll, na.rm = TRUE), mllbound)
# Update profile plot
g1 <- g1 + geom_line(data = data.frame(eta = etagrid, cll = ll-mll),
mapping = aes(x = eta, y = cll), col = "gray") +
  geom_point(data = data.frame(x = 1/max(dat), y = mllbound - mll),
             aes(x = x, y = y),
             col = "gray")


set.seed(202010)
xdat <- evd::rgpd(n = 1000, shape = -0.5)
mle <- mev::fit.pp(xdat, np = 10)$param
fake_mle <- mev::fit.pp(xdat, np = 10,
                        fpar = list(scale = mle['scale']),
                        start = c(-0.1, 1.2))$estimate
fake_mle2 <- evd::fpot(x = xdat,
                       npp = 100,
                       threshold = 0,
                       model = "pp",
                       start = list(loc = -0.1, shape = 1.2),
                       scale = mle['scale'],
                       method = "BFGS")$estimate
pp.score(par = c(fake_mle2[1],mle['scale'],fake_mle2[2]),
         np = 10,
         u = 0,
         dat = xdat)

shape_v <- seq(-1, 1, by = 0.005)
loc_v <- seq(-2,3, by = 0.005)
ll_s <- matrix(NA, nrow = length(shape_v), ncol = length(loc_v))
for(shape_i in seq_along(shape_v)){
  for(loc_i in seq_along(loc_v)){
    ll <- try(pp.ll(par = c(loc_v[loc_i], mle[2], shape_v[shape_i]),
                dat = xdat,
                u = 0, np = 10))
    if(!inherits(ll, "try-error")){
      ll_s[shape_i, loc_i] <- ll
    }
  }
}
hyperbola <- function(shape){
  ifelse(shape > 0,
         mle['scale']/shape,
         ifelse(shape < 0,
         max(xdat) + mle['scale']/shape, NA))
  }
g2 <- ggplot(data = data.frame(loc = rep(loc_v, each = length(shape_v)),
                               shape = rep(shape_v, length.out = length(ll_s)),
                               ll = c(ll_s - max(ll_s, na.rm = TRUE))),
             aes(y = loc, x = shape, z = ll)) +
  geom_raster(aes(fill = exp(ll/10000))) +
  scale_fill_continuous(high = "grey10", low = "grey90", na.value = "transparent") +
  geom_contour(colour = "white", bins = 10, breaks = seq(0, -1e4, by = -1e3)) +
  theme_classic() +
  scale_y_continuous(expand = c(0,0),
                     limits = range(loc_v),
                     # breaks = -2:3,
                     # labels = paste0("$",-2:3,"$")
                     ) +
  scale_x_continuous(expand = c(0,0),
                     limits = range(shape_v),
                     # breaks = seq(-1,1,by=0.5),
                     # labels = paste0("$",c(-1,-0.5,0,0.5,1),"$")
                     ) +
  theme(legend.position = 'none') +
  geom_function(fun = hyperbola, n = 2001) +
  geom_point(mapping = aes(y = mle[1], x = mle[3]), col = "white") +
  geom_point(mapping = aes(y = mle[1], x = mle[3]), pch = 1, col = "black") +
  labs(y = expression(mu), x = expression(xi))
  # labs(y = '$\\mu$', x = '$\\xi$')


g1t <- g1 + labs(
  x = expression(eta),
  y = "profile log likelihood") +
  scale_y_continuous(limits = c(-7,0),
                     breaks = seq(-8L, 0, by = 2L),
                     labels = paste0("$",seq(-8L, 0, by = 2L),"$")) +
  scale_x_continuous(limits = c(-5,1),
                     breaks = seq(-5L, 1L, by=1L),
                     labels = paste0("$",seq(-5L, 1L, by=1L),"$"))

g2t <- g2 +
  labs(y = expression(mu), 
       x = expression(xi)) +
  scale_y_continuous(expand = c(0,0),
                     limits = range(loc_v),
                     breaks = -2:3,
                     labels = paste0("$",-2:3,"$")
  ) +
  scale_x_continuous(expand = c(0,0),
                     limits = range(shape_v),
                     breaks = seq(-1,1,by=0.5),
                     labels = paste0("$",c(-1,-0.5,0,0.5,1),"$")
  )

pdf("fig1_univariate.pdf", width = 9, height = 4)
g1 + g2
dev.off()
