data{
    int N;
    vector[N] l_size;
    vector[N] age;
    array[N] int ind;
}
parameters{
    real<lower=0> alpha;
    vector[N] z;
    real<lower=0> mu;
    real<lower=0> sigma;
    real<lower=0> sig_s;
}
model{
    alpha ~ normal(0,.5);
    z ~ normal(0,1);
    mu ~ normal(0,.5);
    sigma ~ normal(0,.5);
    sig_s ~ normal(0, .5);

    vector[N] s;
    for(i in 1:N){
        s[i] = alpha + (z[ind[i]]*sigma + mu)*age[i];
    }

    l_size ~ normal(s, sig_s);
}