// =====================================================================
//  cdtmbnma.stan
//  Component dose-response network meta-analysis with a dose-dependent
//  interaction surface.
//
//  Arm-based model. Each arm carries a dose for every component (0 = the
//  component is absent). The linear predictor sums a free per-study
//  reference level, a per-component Emax dose-response, and a pairwise
//  dose-dependent interaction surface, with a study-level random effect on
//  active arms. Continuous (normal, known SE) and binary (binomial-logit)
//  outcomes are supported. The interaction surface is additive (none),
//  bilinear, or saturating (general pharmacodynamic interaction).
//
//  This generalises the validated two-component model in cdt_stage1.stan to
//  C components and a binary likelihood. Randomisation is respected through
//  free per-study reference levels; only relative effects are pooled.
//
//  generated quantities returns the pointwise log-likelihood (for loo/WAIC)
//  and the relative effect against the all-zero reference at a user-supplied
//  dose grid (for prediction at unobserved combinations).
// =====================================================================
data {
  int<lower=1> N;                          // arms (one timepoint per arm)
  int<lower=1> S;                          // studies
  int<lower=1> C;                          // components
  int<lower=0> P;                          // interaction pairs

  array[N] int<lower=1, upper=S> study;    // study index of each arm
  matrix<lower=0>[N, C] D;                  // component doses per arm (0 = absent)
  array[N] int<lower=0, upper=1> is_ref;    // 1 = study-baseline arm (no random effect)

  array[P] int<lower=1, upper=C> pair_a;   // first component of each interaction pair
  array[P] int<lower=1, upper=C> pair_b;   // second component of each interaction pair

  vector<lower=0>[C] dstar;                // dose-normalisation references (bilinear surface)

  int<lower=1, upper=2> outcome;           // 1 = continuous, 2 = binary
  int<lower=0, upper=2> interaction;       // 0 = none, 1 = bilinear, 2 = gpdi (saturating)

  // continuous outcome (used when outcome == 1)
  vector[N] y;                             // arm mean change
  vector<lower=0>[N] se;                   // arm standard error (known)

  // binary outcome (used when outcome == 2)
  array[N] int<lower=0> r;                 // events
  array[N] int<lower=0> nn;                // sample size

  // priors
  real<lower=0> prior_emax_sd;
  real prior_logED50_mean;
  real<lower=0> prior_logED50_sd;
  real<lower=0> prior_int_sd;
  real<lower=0> prior_ref_sd;
  real<lower=0> prior_omega_sd;

  // prediction grid
  int<lower=0> G;
  matrix<lower=0>[G > 0 ? G : 1, C] Dpred;
}

parameters {
  vector[S] m;                                       // per-study reference level
  vector[C] emax;                                    // per-component Emax
  vector[C] logED50;                                 // per-component log ED50
  vector[interaction == 1 ? P : 0] eta;              // bilinear interaction (per pair)
  vector[interaction == 2 ? P : 0] INT;              // saturating asymptote (per pair)
  vector<lower=0>[interaction == 2 ? C : 0] kappa;   // saturating half-dose (per component)
  real<lower=0> omega;                               // random-effect SD
  vector[N] z;                                       // non-centred random effects
}

transformed parameters {
  vector<lower=0>[C] ED50 = exp(logED50);
  vector[N] mu;                                       // linear predictor per arm
  for (i in 1:N) {
    real lp = m[study[i]];
    for (c in 1:C)
      lp += emax[c] * D[i, c] / (ED50[c] + D[i, c]);
    if (interaction == 1) {
      for (p in 1:P)
        lp += eta[p]
              * (D[i, pair_a[p]] / dstar[pair_a[p]])
              * (D[i, pair_b[p]] / dstar[pair_b[p]]);
    } else if (interaction == 2) {
      for (p in 1:P) {
        int a = pair_a[p];
        int b = pair_b[p];
        lp += INT[p]
              * (D[i, a] / (kappa[a] + D[i, a]))
              * (D[i, b] / (kappa[b] + D[i, b]));
      }
    }
    if (is_ref[i] == 0)
      lp += z[i] * omega;
    mu[i] = lp;
  }
}

model {
  m       ~ normal(0, prior_ref_sd);
  emax    ~ normal(0, prior_emax_sd);
  logED50 ~ normal(prior_logED50_mean, prior_logED50_sd);
  if (interaction == 1)
    eta ~ normal(0, prior_int_sd);
  if (interaction == 2) {
    INT   ~ normal(0, prior_int_sd);
    kappa ~ lognormal(prior_logED50_mean, prior_logED50_sd);
  }
  omega ~ normal(0, prior_omega_sd);   // half-normal via <lower=0>
  z     ~ std_normal();

  if (outcome == 1)
    y ~ normal(mu, se);
  else
    r ~ binomial_logit(nn, mu);
}

generated quantities {
  vector[N] log_lik;
  for (i in 1:N) {
    if (outcome == 1)
      log_lik[i] = normal_lpdf(y[i] | mu[i], se[i]);
    else
      log_lik[i] = binomial_logit_lpmf(r[i] | nn[i], mu[i]);
  }

  // relative effect vs the all-zero reference (random effect marginalised at 0)
  vector[G] pred;
  for (g in 1:G) {
    real lp = 0;
    for (c in 1:C)
      lp += emax[c] * Dpred[g, c] / (ED50[c] + Dpred[g, c]);
    if (interaction == 1) {
      for (p in 1:P)
        lp += eta[p]
              * (Dpred[g, pair_a[p]] / dstar[pair_a[p]])
              * (Dpred[g, pair_b[p]] / dstar[pair_b[p]]);
    } else if (interaction == 2) {
      for (p in 1:P) {
        int a = pair_a[p];
        int b = pair_b[p];
        lp += INT[p]
              * (Dpred[g, a] / (kappa[a] + Dpred[g, a]))
              * (Dpred[g, b] / (kappa[b] + Dpred[g, b]));
      }
    }
    pred[g] = lp;
  }
}
