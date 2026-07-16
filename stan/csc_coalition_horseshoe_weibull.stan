// Causal Synergy Curve -- regularized-horseshoe coalition surface, Weibull proportional hazards.
//
// Same structured linear predictor as the binary model:
//   eta = alpha0 + Xb*beta_b + Xc*theta   (baseline block Normal, coalition block horseshoe)
// Weibull PH cumulative hazard:  H(t) = t^shape * exp(eta);  S(t)=exp(-H);  F(t)=1-S.
// This reproduces the applied standardization used in the manuscript:
//   F_S(t) = 1 - exp( -exp(eta) * t^shape ).
data {
  int<lower=1> N;
  int<lower=1> Pb;
  int<lower=1> Pc;
  matrix[N, Pb] Xb;
  matrix[N, Pc] Xc;
  vector<lower=0>[N] time;                    // follow-up time (>0)
  array[N] int<lower=0, upper=1> event;       // 1 = event, 0 = right-censored
  real<lower=0> scale_global;
  real<lower=0> slab_scale;
  real<lower=0> slab_df;
  real<lower=0> beta_scale;
}
transformed data {
  vector[N] log_time;
  vector[N] event_v;
  for (i in 1:N) {
    log_time[i] = log(time[i]);
    event_v[i] = event[i];
  }
}
parameters {
  real alpha0;
  vector[Pb] beta_b;
  vector[Pc] z;
  vector<lower=0>[Pc] lambda;
  real<lower=0> tau;
  real<lower=0> caux;
  real<lower=0> shape;                         // Weibull shape kappa
}
transformed parameters {
  vector[Pc] theta;
  {
    real c2 = square(slab_scale) * caux;
    vector[Pc] lam2 = square(lambda);
    vector[Pc] lambda_tilde = sqrt( c2 * lam2 ./ (c2 + square(tau) * lam2) );
    theta = z .* lambda_tilde * tau;
  }
}
model {
  alpha0 ~ normal(0, 5);
  beta_b ~ normal(0, beta_scale);
  z ~ std_normal();
  lambda ~ cauchy(0, 1);
  tau ~ cauchy(0, scale_global);
  caux ~ inv_gamma(0.5 * slab_df, 0.5 * slab_df);
  shape ~ gamma(2, 2);                          // weakly-informative around 1
  {
    vector[N] eta = alpha0 + Xb * beta_b + Xc * theta;
    vector[N] logH = eta + shape * log_time;
    vector[N] logH_safe;
    for (i in 1:N) {
      logH_safe[i] = fmin(logH[i], 30);                    // stable in proposal tails
    }
    target += -sum(exp(logH_safe));                        // log S(t)
    target += dot_product(event_v, rep_vector(log(shape), N) +
                          (shape - 1) * log_time + eta);   // event log h(t)
  }
}
generated quantities {
  vector[Pc] kappa_incl;
  {
    real c2 = square(slab_scale) * caux;
    vector[Pc] lam2 = square(lambda);
    vector[Pc] s2 = square(tau) * (c2 * lam2 ./ (c2 + square(tau) * lam2));
    for (j in 1:Pc) kappa_incl[j] = s2[j] / (s2[j] + 1);
  }
}
