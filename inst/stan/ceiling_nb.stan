// Negative-binomial model of two unit-by-period count records, so that the
// concentration of the EXPECTED series, the object the ceiling actually
// bounds, can be given a posterior. Each record has a unit effect and a
// period effect on the log scale, both partially pooled, and its own
// clumping. The two records share nothing: no lag is fitted here, and the
// expected series are handed back to R, where the ceiling is taken over the
// same kernels as everywhere else in the package.
//
// Non-centred throughout. Priors: the intercepts are centred on the log of
// each record's mean count with a scale of 2.5; the effect scales are
// half-Student-t(3, 0, 1); the clumping is parameterised as 1/sqrt(phi) with
// a half-normal(0, 1) prior, which puts most mass on moderate to strong
// clumping and lets the data pull toward Poisson.
data {
  int<lower=1> H;                    // units
  int<lower=1> T;                    // periods of the reproductive record
  int<lower=1> J;                    // scored recruit periods
  array[H, T] int<lower=0> X;        // reproduction, units by periods
  array[H, J] int<lower=0> R;        // recruits, units by scored periods
  real mx;                           // prior centre of the reproduction intercept
  real mr;                           // prior centre of the recruit intercept
}
parameters {
  real a_x;
  real a_r;
  vector[H] zu_x;
  vector[T] zv_x;
  vector[H] zu_r;
  vector[J] zv_r;
  real<lower=0> s_ux;
  real<lower=0> s_vx;
  real<lower=0> s_ur;
  real<lower=0> s_vr;
  real<lower=0> ir_x;                // 1 / sqrt(phi_x)
  real<lower=0> ir_r;                // 1 / sqrt(phi_r)
}
transformed parameters {
  vector[H] u_x = s_ux * zu_x;
  vector[T] v_x = s_vx * zv_x;
  vector[H] u_r = s_ur * zu_r;
  vector[J] v_r = s_vr * zv_r;
  real phi_x = inv_square(ir_x);
  real phi_r = inv_square(ir_r);
}
model {
  a_x ~ normal(mx, 2.5);
  a_r ~ normal(mr, 2.5);
  zu_x ~ std_normal();
  zv_x ~ std_normal();
  zu_r ~ std_normal();
  zv_r ~ std_normal();
  s_ux ~ student_t(3, 0, 1);
  s_vx ~ student_t(3, 0, 1);
  s_ur ~ student_t(3, 0, 1);
  s_vr ~ student_t(3, 0, 1);
  ir_x ~ normal(0, 1);
  ir_r ~ normal(0, 1);
  for (h in 1:H) {
    X[h] ~ neg_binomial_2_log(a_x + u_x[h] + v_x, phi_x);
    R[h] ~ neg_binomial_2_log(a_r + u_r[h] + v_r, phi_r);
  }
}
generated quantities {
  // expected system totals at each period: the series the ceiling is
  // computed from (EX) and compared with (ER)
  vector[T] EX;
  vector[J] ER;
  for (t in 1:T) EX[t] = sum(exp(a_x + u_x + v_x[t]));
  for (j in 1:J) ER[j] = sum(exp(a_r + u_r + v_r[j]));
}
