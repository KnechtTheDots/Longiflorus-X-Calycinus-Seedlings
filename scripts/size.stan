data{

    // main data

    int N; // total number of observations
    int N_obs; // number that were observed for the 1 Day after germination variable
    int N_miss; // number that were not observed for the 1 day after germination variable
    array[N] int ind; // index variable for individuals
    array[N_obs] int ii_obs; // index variable for the 1dag observed individuals
    array[N_miss] int ii_miss; // index variable for the 1dag not observed individuals
    vector[N_obs] start_obs; // size 1dag for the individuals that were observed
    vector[N] seed_size; // seed size
    vector[N] germ_day; // germination day
    vector[N] size_end; // final size before start of drought treatment
    
    // rgr estimation model

    int N_rgr; // observations in the relative growth rate dataset (more than N because each individual was observed multiple times)
    vector[N_rgr] size_rgr; // log size at each timepoint
    vector[N_rgr] age; // age at each timepoint
    array[N_rgr] int ind_rgr; // index that links to ind above

    //survival modle

    array[N] int survive; // whether or not the individual was alive (1) or dead (0) on day 3 of the drought treatment

    // fake data to predict over

    int N_tilde; // number of data points to predict over
    vector[N_tilde] x_tilde; // dependent varianbles (-3,3) to predict over
}
transformed data{
  vector[N_rgr] rgr_size = log(size_rgr);
}
parameters{

    //rgr model parameters

    vector<lower=0>[N_rgr] alpha_rgr; // intercept for rgr model
    vector<lower=0>[N] rgr; // predicted relative growth rates
    real<lower=0> mu_rgr; // hyperparameter for rgr
    real<lower=0> sigma_rgr; // hyperparameter for rgr
    real<lower=0> sigma_size; // variance for the rgr model
    real<lower=0> sig_a; // hyperparameter for alpha_rgr

    //size model parameters

    vector[N_miss] start_miss; // values of the missing observations
    real alpha; // intercept for the final size model
    real beta_start; // coefficient for the size 1dag
    real beta_seed; // coefficient for the seed size
    real beta_germ; // coefficient for the germination timing
    real beta_rgr; // coefficient for relative growth rate
    real<lower=0> sigma; // standard deviation for the size model
    real bsg;
    real bgs;

    // survival model

    real alpha_surv; // intercept for the survival model
    real beta_surv; // coefficient for the relationship between size and survival
}
transformed parameters{

    // make a vector containing the observed 1 day after germ sizes and the predicted/missing
    vector[N] size_start; 
    size_start[ii_obs] = start_obs;
    size_start[ii_miss] = start_miss;

    // standardize rgr for the final size model
    vector[N] z = (rgr - mean(rgr))/sd(rgr);
}
model{

    //priors for rgr model
    alpha_rgr ~ normal(0, sig_a);
    rgr ~ normal(mu_rgr, sigma_rgr);
    log(mu_rgr) ~ normal(0, 1);
    log(sigma_rgr) ~ normal(0,1);
    log(sigma_size) ~ normal(0, 1);
    log(sig_a) ~ normal(0, 1);

    // RGR model
    vector[N_rgr] s;
    vector[N_rgr] a;
    vector[N_rgr] b;
    for(i in 1:N_rgr){
        s[i] = alpha_rgr[ind_rgr[i]] + rgr[ind_rgr[i]] * age[i];
        //a[i] = (s[i]^2)/sigma_size;
        //b[i] = s[i]/sigma_size;
    }
    
    // [log size | rgr]
    rgr_size ~ normal(s, sigma_size);

    // priors for size model
    start_miss ~ normal(0,1);
    alpha ~ normal(0,.5);
    beta_start ~ normal(0,.5);
    beta_seed ~ normal(0, .5);
    beta_germ ~ normal(0, .5);
    beta_rgr ~ normal(0, .5);
    log(sigma) ~ normal(0, 1);
    bsg ~ normal(0, .5);
    bgs ~ normal(0, .5);
    vector[N] mu;
    for(i in 1:N){
        mu[i] = alpha + (beta_start + bsg*germ_day[i])*size_start[i] + (beta_seed + bgs*germ_day[i])*seed_size[i] + beta_germ*germ_day[i] + beta_rgr*z[ind[i]];
    }

    // [final size | rgr, seed size, germ dat, size 1 day after germ]
    size_end ~ normal(mu, sigma);

    // survival model
    alpha_surv ~ normal(0, .5);
    beta_surv ~ normal(0, .5);

    // [survive | size]
    survive ~ bernoulli_logit(alpha_surv + beta_surv*mu);
}
generated quantities {
    vector[N] log_lik;
    for(i in 1:N){
        log_lik[i] = bernoulli_logit_lpmf(survive[i]|alpha_surv + beta_surv*size_end[i]);
    }

    // predict the final size of a plant for posterior predictive checks
    array[N] real y_rep = normal_rng(alpha + (beta_start + bsg*germ_day) .* size_start + (beta_seed + bgs*germ_day).*seed_size + beta_germ*germ_day + beta_rgr*z, sigma);
    

    // posterior of prob of survival over size (-3,3)
    vector[N_tilde] p_tilde = inv_logit(alpha_surv + beta_surv * x_tilde);

    // posterior predictive of rgr
    array[N_tilde] real y_t_rgr = normal_rng(alpha + beta_rgr * x_tilde, sigma); // predict size with everything else held at the mean

    vector[N_tilde] p_rgr = inv_logit(alpha_surv + beta_surv * to_vector(y_t_rgr));
}