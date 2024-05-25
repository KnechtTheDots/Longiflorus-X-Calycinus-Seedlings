bayes_boot <- function(x){
  weights_mat <- gtools::rdirichlet(n = 1e4, alpha = rep(1, length(x)))
  mu <- apply(weights_mat, 1, function(z) sum(z*x))
  sigma <- apply(weights_mat, 1, function(z) sqrt(sum(z*(x - mean(x))^2)))
  return(data.frame(mu = mu, sd = sigma))
}


post_pred <- function(df, trait){
  d <- df %>% 
    filter(germ_day==3) %>% 
    drop_na(trait) 
  
  lon <- d %>% 
    filter(class == "LON") %>% 
    pull(trait)
  
  cal <- d %>% 
    filter(class == "CAL") %>% 
    pull(trait)
  
  f1 <- d %>% 
    filter(class == "F1") %>% 
    pull(trait)
  
  f2 <- d %>% 
    filter(class == "F2") %>% 
    pull(trait)
  
  mid_p <- (mean(lon)+mean(cal))/2
  
  lon <- lon/mid_p
  cal <- cal/mid_p
  f1 <- f1/mid_p
  f2 <- f2/mid_p
  
  library(gtools)
  
  lon_boot <- bayes_boot(lon)
  cal_boot <- bayes_boot(cal)
  f1_boot <- bayes_boot(f1)
  f2_boot <- bayes_boot(f2)
  
  lon_pp <- matrix(nrow=1e4, ncol = length(lon))
  cal_pp <- matrix(nrow=1e4, ncol = length(cal))
  f1_pp <- matrix(nrow=1e4, ncol = length(f1))
  f2_pp <- matrix(nrow=1e4, ncol = length(f2))
  
  delta_pp <- c()
  for(i in 1:1e4){
    lon_pp[i,] <- rnorm(length(lon), lon_boot$mu[i], lon_boot$sd[i])
    cal_pp[i,] <- rnorm(length(cal), cal_boot$mu[i], cal_boot$sd[i])
    f1_pp[i,] <- rnorm(length(f1), f1_boot$mu[i], f1_boot$sd[i])
    f2_pp[i,] <- rnorm(length(f2), f2_boot$mu[i], f2_boot$sd[i])
    
    delta_pp[i] <- mean(f2_pp[i,]) - (.25*(mean(lon_pp[i,]) + mean(cal_pp[i,])) + .5*mean(f1_pp[i,]))
  }
  return(list(lon_pp = lon_pp, cal_pp = cal_pp, f1_pp = f1_pp, f2_pp = f2_pp,
              delta_pp = delta_pp))
}