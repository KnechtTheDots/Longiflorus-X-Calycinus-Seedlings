data{
    int N;
    vector[N] l_size;
    vector[N] age;
    array[N] int ind;
}
parameters{
    vector<lower=0>[N] alpha_rgr;
    vector<lower=0>[N] rgr;
    real<lower=0> mu_rgr;
    real<lower=0> sigma_rgr;
    real<lower=0> sigma_size;
    //real<lower=0> mu_a;
    real<lower=0> sig_a;
    //vector[N] z_a;
}
model{
    rgr ~ normal(mu_rgr, sigma_rgr);
    alpha_rgr ~ normal(.25, sig_a);
    mu_rgr ~ exponential(3);
    sigma_rgr ~ exponential(3);
    sigma_size ~ normal(.5, 1);
    //mu_a ~ exponential(3);
    sig_a ~ exponential(3);
    

    vector[N] s;
    vector[N] a;
    vector[N] b;
    for(i in 1:N){
        s[i] = alpha_rgr[ind[i]] + rgr[ind[i]] * age[i];
        a[i] = (s[i]^2)/sigma_size;
        b[i] = s[i]/sigma_size;
    }

    l_size ~ gamma(a,b);
}
generated quantities {
    vector[N] s;
    vector[N] a;
    vector[N] b;
    for(i in 1:N){
        s[i] = alpha_rgr[ind[i]] + rgr[ind[i]] * age[i];
        a[i] = (s[i]^2)/sigma_size;
        b[i] = s[i]/sigma_size;
    }
    array[N] real size_rep = gamma_rng(a,b);
}