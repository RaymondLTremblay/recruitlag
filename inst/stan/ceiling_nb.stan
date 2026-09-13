// Negative-binomial model of the recruit counts on every unit at every scored
// period, with a unit effect and a period effect on the log scale, so that
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
data {
  int<lower=1> H;                    // units
  int<lower=1> J;                    // scored recruit periods
  array[H, J] int<lower=0> R;        // recruits, units by scored periods
  real mr;                           // prior centre of the intercept (log scale)
  real<lower=0> sigma_scale;         // scale of the half-t priors on the effect sds
  real<lower=0> phi_shape;           // gamma prior on the clumping parameter
  real<lower=0> phi_rate;
}
transformed data {
  // the counts as one flat vector, with the unit and period of each cell, so
  // that the likelihood is a single vectorised call rather than a loop
  int N = H * J;
  array[N] int<lower=0> Rf;
  array[N] int<lower=1, upper=H> hh;
  array[N] int<lower=1, upper=J> jj;
  for (h in 1:H) for (j in 1:J) {
    int n = (h - 1) * J + j;
    Rf[n] = R[h, j]; hh[n] = h; jj[n] = j;
  }
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
  Rf ~ neg_binomial_2_log(a_r + u_r[hh] + v_r[jj], phi_r);
}
generated quantities {
  // expected system total at each scored period: the series whose
  // concentration is compared with the ceiling
  vector[J] ER;
  for (j in 1:J) ER[j] = sum(exp(a_r + u_r + v_r[j]));
}
