data {
    int N;
    int N_obs;
    int N_miss;
    vector[N] seed_size;
    vector[N_obs] size_start;
    array[N_obs] int ii_obs;
    array[N_miss] int ii_miss;
    vector[N] age;
    vector[N] size_end;
}
parameters{
    // size_start model
    vector[N_miss] size_miss;
    real alpha_start;
    real beta_seed;
    real<lower=0> sigma_start;

    // final size model
    //vector<lower=0>[N] rgr;
    //vector[N] z;
    //real<lower=0> mu_rgr;
    //real sigma_rgr;
    real alpha_end;
    real beta_start;
    real<lower=0> sigma_end;
}
transformed parameters {
   vector[N] size_0;
   size_0[ii_obs] = size_start;
   size_0[ii_miss] = size_miss;
}
model{
    // size_start model
    alpha_start ~ normal(0, 1);
    beta_seed ~ normal(0, 1);
    sigma_start ~ exponential(1);
    vector[N] mu_start;
    for(i in 1:N){
        mu_start[i] = alpha_start + beta_seed * seed_size[i];
    }

    size_start ~ normal(mu_start[ii_obs], sigma_start);
    size_miss ~ normal(mu_start[ii_miss], sigma_start);

    //final size model
        // rgr
    //rgr ~ normal(0, 1);
    //z ~ normal(0, 1);
    //mu_rgr ~ normal(.25, .5);
    //sigma_rgr ~ normal(0, 1);
    alpha_end ~ normal(0, 1);
    beta_start ~ normal(0, 1);
    sigma_end ~ exponential(1);

    vector[N] mu_end;
    for(i in 1:N){
        mu_end[i] = alpha_end + beta_start * size_0[i];// + rgr[i] * age[i];
    }

    size_end ~ normal(mu_end, sigma_end);
}