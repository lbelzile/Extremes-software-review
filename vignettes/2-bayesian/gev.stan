functions {
  real gev_lpdf(vector y, vector mu, vector sigma, real xi) {
   // generalized extreme value log density
   real xis = xi;
   if(abs(xi) < 1e-3){
    xis = (2*step(xi) - 1)*1e-3;
   }
   real inv_xi = inv(xis);
   int N = rows(y);
   vector[N] z = (y - mu) ./ sigma;
      if(xi < 0 && max(z) >= -inv_xi)
      	reject("Outside of bounds with xi negative", xi);
      if(xi > 0 && min(z) <= -inv_xi)
      reject("Outside of bounds with xi positive", xi);
    if (min(sigma) <= 0)
      reject("sigma<=0; found sigma =", min(sigma));
    if (abs(xi) > 1e-3)
      return -(1+inv_xi)*sum(log1p(xi*z)) -sum(log(sigma)) - sum(exp(-inv_xi*log1p(xi*z)));
    else
     return -sum(log(sigma)) + (0.3 - abs(xi))/0.3 * (-(1+inv_xi)*sum(log1p(xis*z)) - sum(exp(-inv_xi*log1p(xis*z)))) +
     (abs(xi))/0.3 * (-sum(z) - sum(exp(-z))); // linear interpolation around xi=0
  }
  
  real gev_rng(real mu, real sigma, real xi) {
    // generalized extreme value random number generation
    if (sigma<=0)
      reject("sigma<=0; found sigma =", sigma);
    if (abs(xi) > 1e-15)
      return mu + sigma*(exponential_rng(1)^xi-1)/xi;
    else
      return mu - sigma*log(exponential_rng(1));  // xi=0
  }
}

data {
  int<lower=0> N;
  int<lower=0> p_loc;
  int<lower=0> p_scale;
  vector[N] y;
  vector[3] priorvar;
  matrix[N, p_loc] X_loc;
  matrix[N, p_scale] X_scale;
}
parameters {
  vector[p_loc] beta_loc;
  vector[p_scale] beta_scale;
  real<lower=-1> xi; // maximum likelihood is infinite beyond
}
model {
  xi ~ normal(0, priorvar[3]);
  beta_loc ~ normal(0, priorvar[1]);
  beta_scale ~ normal(0, priorvar[2]);
  y ~ gev(X_loc * beta_loc, exp(X_scale * beta_scale), xi);
}

