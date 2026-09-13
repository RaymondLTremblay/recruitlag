// Negative-binomial model of the recruit counts on every unit at every scored
// period at which that unit was censused, with a unit effect and a period effect on the log scale, so that
// the concentration of the EXPECTED recruitment series, the object the
// ceiling actually bounds, can be given a posterior. Only the recruits are
// modelled: the reproductive record is the driver the hypothesis names, and
// it enters the ceiling as observed. (An earlier version modelled both
// records; when the reproductive period effects shrank toward zero the
// modelled ceiling went to zero and the exceedance, a ratio, exploded.)
//
// Non-centred throughout. The priors are the working priors of the companion
// paper's hierarchical recruitment model (phase 4d): intercept normal(centre,
// 2), effect scales half-Student-t(3, 0, sigma_scale), and the clumping
// phi ~ gamma(phi_shape, phi_rate), gamma(2, 0.1) by default. They are
// placeholders pending elicitation, and the three hyperparameters are data so
// that a user can change them without editing the model.
functions {
  // one slice of the likelihood, for reduce_sum (within-chain threads)
  real partial_nb(array[] int Rs, int start, int end, vector eta, real phi) {
    return neg_binomial_2_log_lpmf(Rs | eta[start:end], phi);
  }
}
data {
  int<lower=1> H;                    // units
  int<lower=1> J;                    // scored recruit periods
  int<lower=1> N;                    // observed unit-by-period cells (a cell not censused is absent, not zero)
  array[N] int<lower=0> Rf;          // recruit count in each observed cell
  array[N] int<lower=1, upper=H> hh; // its unit
  array[N] int<lower=1, upper=J> jj; // its period
  real mr;                           // prior centre of the intercept (log scale)
  real<lower=0> sigma_scale;         // scale of the half-t priors on the effect sds
  real<lower=0> phi_shape;           // gamma prior on the clumping parameter
  real<lower=0> phi_rate;
  int<lower=0> grainsize;            // 0: one vectorised call; > 0: reduce_sum slices of about this size
}
parameters {
  real a_r;
  vector[H] zu_r;
  vector[J] zv_r;
  real<lower=0> s_ur;
  real<lower=0> s_vr;
  real<lower=0> phi_r;
}
transformed parameters {
  vector[H] u_r = s_ur * zu_r;
  vector[J] v_r = s_vr * zv_r;
}
model {
  a_r ~ normal(mr, 2);
  zu_r ~ std_normal();
  zv_r ~ std_normal();
  s_ur ~ student_t(3, 0, sigma_scale);
  s_vr ~ student_t(3, 0, sigma_scale);
  phi_r ~ gamma(phi_shape, phi_rate);
  {
    vector[N] eta = a_r + u_r[hh] + v_r[jj];
    if (grainsize > 0) target += reduce_sum(partial_nb, Rf, grainsize, eta, phi_r);
    else Rf ~ neg_binomial_2_log(eta, phi_r);
  }
}
generated quantities {
  // expected system total at each scored period: the series whose
  // concentration is compared with the ceiling
  vector[J] ER;
  for (j in 1:J) ER[j] = sum(exp(a_r + u_r + v_r[j]));
}
