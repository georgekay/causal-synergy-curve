// Causal Synergy Curve -- regularized-horseshoe coalition surface, binary outcome (cloglog link).
//
// Estimating model (locked 2026-07-09):
//   g{E(Y|W,D)} = alpha0 + h(W) + sum_j beta_j D_j + f(m) + sum_g kappa_g Gate_g + sum_T theta_T phi_T
// The baseline block Xb = [ W, main effects D_j, burden-depth (m, m^2), anchored gates Gate_g ]
// receives a weakly-informative Normal prior. The coalition block Xc = [ products phi_T for
// |T|>=2 up to the chosen order ] receives a REGULARIZED HORSESHOE (Piironen & Vehtari 2017),
// so null coalitions are crushed while genuine ones escape -- the estimand-first, sparsity-aware
// prior the Synergy Curve requires. Link is cloglog to match a proportional-hazard-at-horizon DGP.
data {
  int<lower=1> N;
  int<lower=1> Pb;                 // baseline columns (covariates + mains + depth + gates)
  int<lower=1> Pc;                 // coalition columns (products)
  matrix[N, Pb] Xb;
  matrix[N, Pc] Xc;
  array[N] int<lower=0, upper=1> y;
  real<lower=0> scale_global;      // tau0: global shrink scale, e.g. p0/((Pc-p0))*(1/sqrt(N)) with p0 = expected # active
  real<lower=0> slab_scale;        // c: slab scale (e.g. 2)
  real<lower=0> slab_df;           // slab dof (e.g. 4)
  real<lower=0> beta_scale;        // baseline weakly-informative sd (e.g. 2.5)
}
parameters {
  real alpha0;
  vector[Pb] beta_b;
  vector[Pc] z;
  vector<lower=0>[Pc] lambda;      // local scales
  real<lower=0> tau;               // global scale
  real<lower=0> caux;              // slab auxiliary
}
transformed parameters {
  vector[Pc] theta;                // coalition coefficients (regularized horseshoe)
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
  {
    vector[N] eta = alpha0 + Xb * beta_b + Xc * theta;
    for (i in 1:N) {
      // Stable cloglog Bernoulli likelihood. The upper clamp only affects
      // impossible proposal tails, where exp(eta) is already effectively infinite.
      real H = exp(fmin(eta[i], 30));
      if (y[i] == 1) target += log1m_exp(-H);
      else target += -H;
    }
  }
}
generated quantities {
  vector[Pc] kappa_incl;   // per-coalition posterior inclusion weight in [0,1] (~1 kept, ~0 shrunk)
  {
    real c2 = square(slab_scale) * caux;
    vector[Pc] lam2 = square(lambda);
    vector[Pc] s2 = square(tau) * (c2 * lam2 ./ (c2 + square(tau) * lam2));
    for (j in 1:Pc) kappa_incl[j] = s2[j] / (s2[j] + 1);
  }
}
