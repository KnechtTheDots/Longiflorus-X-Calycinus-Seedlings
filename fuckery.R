d <- read.csv("data/seed_growth_vs_seed_size.csv")
library(tidyverse)


rgr_dat <- d %>% 
  drop_na(germ_day, death_day) %>% 
  mutate(ind = 1:105) %>% 
  pivot_longer(4:10, names_to = "day", values_to = "size") %>% 
  drop_na(size) %>% 
  mutate(l_size = log(size)) %>% 
  separate(day, into = (c("a", "day"))) %>% 
  mutate(day = as.numeric(day),
         age = day - germ_day) %>% 
  filter(age > 0) %>% 
  select(ind, age, l_size)

dat <- list(
  N = nrow(rgr_dat),
  l_size = rgr_dat$l_size,
  ind = rgr_dat$ind,
  age = rgr_dat$age
)

library(cmdstanr)
mod <- cmdstan_model("scripts/rgr_est.stan", compile = F)
mod$check_syntax()
mod <- cmdstan_model("scripts/rgr_est.stan")

fit <- mod$sample(
  data = dat,
  chains = 4,
  parallel_chains = 4
)

y <- rgr_dat$l_size
y_rep <- fit$draws("size_rep", format = "df")[,1:length(y)]


library(bayesplot)

post <- as.array(fit$draws())
dimnames(post)

plot(apply(y_rep, 2, mean)~y)
abline(a = 0, b = 1)





d <- read.csv("data/seed_growth_vs_seed_size.csv")


rgr_dat <- d %>% 
  drop_na(germ_day, death_day) %>% 
  mutate(ind = 1:105) %>% 
  pivot_longer(4:10, names_to = "day", values_to = "size") %>% 
  drop_na(size) %>% 
  mutate(l_size = log(size)) %>% 
  separate(day, into = (c("a", "day"))) %>% 
  mutate(day = as.numeric(day),
         age = day - germ_day) %>% 
  filter(age > 0) %>% 
  select(ind, age, l_size)

full <- d %>% drop_na(germ_day, death_day) %>%
  mutate(ind = 1:105) %>% 
  pivot_longer(4:10, names_to = "day", values_to = "size") %>% 
  drop_na(size) %>%
  drop_na(size) %>% 
  separate(day, c("fubar", "day")) %>% 
  mutate(day = as.numeric(day)) %>% 
  select(ind, seed_size, germ_day, death_day, day, size) %>% 
  mutate(age = day - germ_day) %>%
  filter(age == 1 | day == 17) %>%
  mutate(tp = ifelse(age == 1, "start","finish")) %>% 
  pivot_wider(id_cols = c(ind, seed_size, germ_day, death_day),
              names_from = tp, values_from = size)

mu <- mean(full$start[!(is.na(full$start))])

sd <- sd(full$start[!(is.na(full$start))])

full$start <- ifelse(is.na(full$start), -100, (full$start - mu)/sd)
full$finish <- stn(full$finish)
full$seed_size <- stn(full$seed_size)
full$survive <- full$death_day - 3



stan_dat <- list(
  N = nrow(full),
  N_tilde = 200,
  N_obs = sum(full$start != -100),
  N_miss = sum(full$start == -100),
  ind = full$ind,
  ii_obs = full$ind[full$start != -100],
  ii_miss = full$ind[full$start == -100],
  start_obs = full$start[full$start!=-100],
  size_end = full$finish,
  seed_size = full$seed_size,
  germ_day = stn(full$germ_day),
  N_rgr = nrow(rgr_dat),
  l_size = rgr_dat$l_size,
  age = rgr_dat$age,
  ind_rgr = rgr_dat$ind,
  survive = full$survive,
  x_tilde = seq(-3,3, l = 200)
)




mod <- cmdstan_model("scripts/size.stan", compile = F)
mod$check_syntax()
mod <- cmdstan_model("scripts/size.stan")

fit <- mod$sample(
  data = stan_dat,
  chains = 4, 
  parallel_chains = 4
)


y <- stan_dat$size_end
y_rep <- fit$draws("y_rep", format = "df")[,1:length(y)]

mu <- apply(y_rep, 2, mean)

plot(mu ~ y)
abline(a = 0, b = 1)

post <- as.array(fit$draws())

mcmc_intervals(post, pars = c("bsg", "bgs", "beta_seed", "beta_rgr",
                              "beta_germ", "beta_start", "beta_surv"))

df <- full %>% 
  mutate(age = 17 - germ_day,
         seed_size = log(seed_size),
         start = ifelse(is.na(start), -100, log(start)),
         finish = log(finish))

dat <- list(
  N = nrow(df),
  N_obs = sum(df$start!=-100),
  N_miss = sum(df$start==-100),
  seed_size = df$seed_size,
  size_start = df$start[df$start!=-100],
  ii_obs = df$ind[df$start!=-100],
  ii_miss = df$ind[df$start==-100],
  age = df$age,
  size_end = df$finish
)

mod <- cmdstan_model("scripts/luxury.stan")

fit <- mod$sample(
  data = dat,
  chains = 4,
  parallel_chains = 4
)

post <- as.array(fit$draws(c("z", "mu_rgr", "sigma_rgr", "sigma_end")))

library(bayesplot)
mcmc_parcoord(post)
library(shinystan)
launch_shinystan(fit)

mcmc_pairs(post, pars = c("sigma_rgr", "sigma_end"))
