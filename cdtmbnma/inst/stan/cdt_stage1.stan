// =====================================================================
//  cdtMBNMA -- Stage 1 canonical Stan model
//  Component x Dose x Time Model-Based Network Meta-Analysis
//
//  Two components (A, B); exponential time-course; bilinear (one-parameter)
//  dose-dependent interaction on the asymptote; additive log-linear dose
//  effects on the onset rate; independent between-study random effects on
//  relative effects; univariate-normal aggregate likelihood (within-arm
//  correlation across time = 0).
//
//  Mirrors cdtmbnma/model_numpyro.py exactly. Designed to fork alongside
//  the Strathe et al. (2026) DT-MBNMA Stan code. Randomisation is respected
//  via free per-study reference levels (mE, mk); only relative (component /
//  dose / interaction) effects are pooled.
//
//  generated quantities returns pointwise log-lik (for loo/ELPD) and the
//  long-term combination effect vs placebo at a user-supplied dose grid
//  (for prediction / synergy surfaces / held-out combinations).
// =====================================================================
data {
  int<lower=1> n_obs;                 // rows: study x arm x timepoint
  int<lower=1> n_arms;                // distinct (study, arm) pairs
  int<lower=1> n_studies;             // studies
  int<lower=1> n_active;              // active (non-reference) arms

  // row-level
  array[n_obs] int<lower=1, upper=n_arms> row_arm;   // arm of each row
  vector[n_obs] time;                 // weeks
  vector[n_obs] ybar;                 // observed arm mean
  vector<lower=0>[n_obs] se;          // known standard error

  // arm-level covariates
  array[n_arms] int<lower=1, upper=n_studies> arm_sid;
  vector<lower=0>[n_arms] arm_dA;     // dose of component A (0 = absent)
  vector<lower=0>[n_arms] arm_dB;     // dose of component B (0 = absent)
  array[n_active] int<lower=1, upper=n_arms> active_idx;  // which arms get a RE

  real<lower=0> dAstar;               // dose-normalisation reference, comp A
  real<lower=0> dBstar;               // dose-normalisation reference, comp B

  // prediction grid for generated quantities
  int<lower=0> n_pred;
  vector<lower=0>[n_pred > 0 ? n_pred : 1] pred_dA;
  vector<lower=0>[n_pred > 0 ? n_pred : 1] pred_dB;
}

parameters {
  // study-specific reference (placebo) levels -- free fixed effects
  vector[n_studies] mE;               // reference asymptote
  vector[n_studies] mk;               // reference log onset-rate

  // component dose-response on asymptote (Emax)
  real a1_EmaxA;
  real a2_EmaxB;
  real logED50A;
  real logED50B;

  // dose-dependent interaction (bilinear; 0 = additivity)
  real eta;

  // log-rate dose effects (additive)
  real b1_rateA;
  real b2_rateB;

  // between-study RE SDs (independent in Stage 1)
  real<lower=0> omega_E;
  real<lower=0> omega_k;

  // non-centred REs for active arms
  vector[n_active] zE;
  vector[n_active] zk;
}

transformed parameters {
  real<lower=0> ED50A = exp(logED50A);
  real<lower=0> ED50B = exp(logED50B);

  vector[n_arms] uE = rep_vector(0.0, n_arms);
  vector[n_arms] uk = rep_vector(0.0, n_arms);
  for (a in 1:n_active) {
    uE[active_idx[a]] = zE[a] * omega_E;
    uk[active_idx[a]] = zk[a] * omega_k;
  }

  // per-arm asymptote E and log-rate k
  vector[n_arms] E_arm;
  vector[n_arms] k_arm;
  {
    vector[n_arms] gA = a1_EmaxA * arm_dA ./ (ED50A + arm_dA);
    vector[n_arms] gB = a2_EmaxB * arm_dB ./ (ED50B + arm_dB);
    vector[n_arms] inter = eta * (arm_dA / dAstar) .* (arm_dB / dBstar);
    for (j in 1:n_arms) {
      E_arm[j] = mE[arm_sid[j]] + gA[j] + gB[j] + inter[j] + uE[j];
      k_arm[j] = mk[arm_sid[j]]
                 + b1_rateA * log(arm_dA[j] + 1.0)
                 + b2_rateB * log(arm_dB[j] + 1.0)
                 + uk[j];
    }
  }
}

model {
  // ---- priors (weakly informative) ----
  mE ~ normal(0, 5);
  mk ~ normal(-3, 1);
  a1_EmaxA ~ normal(0, 10);
  a2_EmaxB ~ normal(0, 10);
  logED50A ~ normal(0, 1);
  logED50B ~ normal(0, 1);
  eta ~ normal(0, 5);                 // centred at additivity
  b1_rateA ~ normal(0, 1);
  b2_rateB ~ normal(0, 1);
  omega_E ~ normal(0, 2);             // half-normal via <lower=0>
  omega_k ~ normal(0, 0.5);
  zE ~ std_normal();
  zk ~ std_normal();

  // ---- aggregate likelihood (exponential time-course; SE known) ----
  {
    vector[n_obs] mu;
    for (i in 1:n_obs) {
      int j = row_arm[i];
      mu[i] = E_arm[j] * (1.0 - exp(-exp(k_arm[j]) * time[i]));
    }
    ybar ~ normal(mu, se);
  }
}

generated quantities {
  // pointwise log-likelihood for loo / ELPD model comparison
  vector[n_obs] log_lik;
  for (i in 1:n_obs) {
    int j = row_arm[i];
    real mu = E_arm[j] * (1.0 - exp(-exp(k_arm[j]) * time[i]));
    log_lik[i] = normal_lpdf(ybar[i] | mu, se[i]);
  }

  // long-term (asymptote) combination effect vs placebo on a dose grid
  // RE marginalised at 0; reference level excluded -> pure relative effect
  vector[n_pred] pred_effect;
  for (g in 1:n_pred) {
    real gA = a1_EmaxA * pred_dA[g] / (ED50A + pred_dA[g]);
    real gB = a2_EmaxB * pred_dB[g] / (ED50B + pred_dB[g]);
    real it = eta * (pred_dA[g] / dAstar) * (pred_dB[g] / dBstar);
    pred_effect[g] = gA + gB + it;
  }
}
