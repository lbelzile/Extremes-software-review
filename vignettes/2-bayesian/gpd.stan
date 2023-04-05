functions {
  real gpareto_lpdf(vector y, real thresh, real xi, real sigma) {
    // generalised Pareto log pdf
    int N = rows(y);
    real inv_xi = inv(xi);
    if (xi<0 && max(y-thresh)/sigma > -inv_xi)
      reject("xi<0 and max(y-thresh)/sigma > -1/xi; found xi, sigma =", xi, sigma);
    if (sigma<=0)
      reject("sigma<=0; found sigma =", sigma);
    if (abs(xi) > 1e-15)
      return -(1+inv_xi)*sum(log1p((y-thresh) * (xi/sigma))) -N*log(sigma);
    else
      return -sum(y-thresh)/sigma -N*log(sigma); // limit xi->0
  }
  real gpareto_cdf(vector y, real thresh, real xi, real sigma) {
    // generalised Pareto cdf
    real inv_xi = inv(xi);
    if (xi<0 && max(y-thresh)/sigma > -inv_xi)
      reject("xi<0 and max(y-thresh)/sigma > -1/xi; found xi, sigma =", xi, sigma);
    if (sigma<=0)
      reject("sigma<=0; found sigma =", sigma);
    if (abs(xi) > 1e-15)
      return exp(sum(log1m_exp((-inv_xi)*(log1p((y-thresh) * (xi/sigma))))));
    else
      return exp(sum(log1m_exp(-(y-thresh)/sigma))); // limit xi->0
  }
  real gpareto_lcdf(vector y, real thresh, real xi, real sigma) {
    // generalised Pareto log cdf
    real inv_xi = inv(xi);
    if (xi<0 && max(y-thresh)/sigma > -inv_xi)
      reject("xi<0 and max(y-thresh)/sigma > -1/xi; found xi, sigma =", xi, sigma);
    if (sigma<=0)
      reject("sigma<=0; found sigma =", sigma);
    if (abs(xi) > 1e-15)
      return sum(log1m_exp((-inv_xi)*(log1p((y-thresh) * (xi/sigma)))));
    else
      return sum(log1m_exp(-(y-thresh)/sigma)); // limit xi->0
  }
  real gpareto_lccdf(vector y, real thresh, real xi, real sigma) {
    // generalised Pareto log ccdf
    real inv_xi = inv(xi);
    if (xi<0 && max(y-thresh)/sigma > -inv_xi)
      reject("xi<0 and max(y-thresh)/sigma > -1/xi; found xi, sigma =", xi, sigma);
    if (sigma<=0)
      reject("sigma<=0; found sigma =", sigma);
    if (abs(xi) > 1e-15)
      return (-inv_xi)*sum(log1p((y-thresh) * (xi/sigma)));
    else
      return -sum(y-thresh)/sigma; // limit xi->0
  }
  real gpareto_rng(real thresh, real xi, real sigma) {
    // generalised Pareto rng
    if (sigma<=0)
      reject("sigma<=0; found sigma =", sigma);
    if (abs(xi) > 1e-15)
      return thresh + (uniform_rng(0,1)^-xi -1) * sigma / xi;
    else
      return thresh - sigma*log(uniform_rng(0,1)); // limit xi->0
  }
}
data {
  real thresh;
  int<lower=0> N;
  vector<lower=thresh>[N] y;
}
transformed data {
  real ymax = max(y);
}
parameters {
  real<lower=0> sigma;
  real<lower=-sigma/(ymax-thresh)> xi;
}
model {
  y ~ gpareto(thresh, xi, sigma);
}
