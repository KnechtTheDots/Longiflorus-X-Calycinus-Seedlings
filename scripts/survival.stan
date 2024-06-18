data{
  int N; // number of observations
  array[N] int y; // survived/died on day 3
  matrix[N,2] X; // observations of size
  int K;
  matrix[K,2] X_rep; // for posterior predictive
}
parameters{
  vector[2] beta;
}
model{
  beta ~ normal(0, .5);
  y ~ bernoulli_logit(X*beta);
}
generated quantities{
  vector[K] p = inv_logit(X_rep*beta);
}