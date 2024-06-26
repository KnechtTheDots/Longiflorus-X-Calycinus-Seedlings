d <- read.csv("data/seed_growth_vs_seed_size.csv")

d %>% drop_na(day_17, death_day) %>% 
  mutate(size = (day_17 - mean(day_17))/sd(day_17)) %>% 
  mutate(survive = death_day - 3)

dat <- d %>% 
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

stan_dat <- list(
  N = nrow(dat),
  age = dat$age,
  ind = dat$ind,
  l_size = dat$l_size
)

library(cmdstanr)
mod <- cmdstan_model("scripts/rgr.stan")


fit <- mod$sample(
  dat = stan_dat,
  chains = 4,
  parallel_chains = 4
)


library(tidyverse)
d <- read.csv("data/seed_growth_vs_seed_size.csv")

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


X <- matrix(c(rep(1, nrow(full)),
              full$seed_size,
              full$germ_day), ncol = 3)

stan_dat <- list(
  N = nrow(full),
  N_obs = sum(full$start != -100),
  N_miss = sum(full$start == -100),
  ind = full$ind,
  ii_obs = full$ind[full$start != -100],
  ii_miss = full$ind[full$start == -100],
  start_obs = full$start[full$start!=-100],
  size_end = full$finish,
  seed_size = full$seed_size,
  germ_day = full$germ_day,
  N_rgr = nrow(rgr_dat),
  l_size = rgr_dat$l_size,
  age = rgr_dat$age,
  ind_rgr = rgr_dat$ind,
  survive = full$survive
)


library(cmdstanr)
mod <- cmdstan_model("scripts/size.stan", compile = F)
mod$check_syntax()
mod <- cmdstan_model("scripts/size.stan")


fit <- mod$sample(
  data = stan_dat,
  chains = 4,
  parallel_chains = 4
)

fit$print(max_rows = 336)



y <- stan_dat$size_end
library(shinystan)
launch_shinystan(fit)

