# z score function
stn <- function(x){
  (x - mean(x))/sd(x)
}



## functionto do the bayesian bootstrap and get mean and sd
bayes_boot <- function(x){
  weights_mat <- gtools::rdirichlet(n = 1e4, alpha = rep(1, length(x)))
  mu <- apply(weights_mat, 1, function(z) sum(z*x))
  sigma <- apply(weights_mat, 1, function(z) sqrt(sum(z*(x - sum(z*x))^2)))
  return(data.frame(mu = mu, sd = sigma))
}


## generate the posterior predictive distributions

post_pred <- function(df, trait, lik){
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
  
  # mid_p <- (mean(lon)+mean(cal))/2
  # 
  # lon <- lon/mid_p
  # cal <- cal/mid_p
  # f1 <- f1/mid_p
  # f2 <- f2/mid_p
  
  mu <- mean(f2)
  s <- sd(f2)
  
  lon <- (lon - mu)/s
  cal <- (cal - mu)/s
  f1 <- (f1 - mu)/s
  f2 <- (f2 - mu)/s
  
  library(gtools)
  
  lon_boot <- bayes_boot(lon)
  cal_boot <- bayes_boot(cal)
  f1_boot <- bayes_boot(f1)
  f2_boot <- bayes_boot(f2)
  e_f2_boot <- .25*(lon_boot$mu + cal_boot$mu) + .5*f1_boot$mu
  
  lon_pp <- matrix(nrow=1e4, ncol = length(lon))
  cal_pp <- matrix(nrow=1e4, ncol = length(cal))
  f1_pp <- matrix(nrow=1e4, ncol = length(f1))
  f2_pp <- matrix(nrow=1e4, ncol = length(f2))
  
  delta_pp <- c()
  for(i in 1:1e4){
    if(lik=="normal"){
      lon_pp[i,] <- rnorm(length(lon), mean = lon_boot$mu[i], 
                          sd = lon_boot$sd[i])
      cal_pp[i,] <- rnorm(length(cal), mean = cal_boot$mu[i], 
                          sd = cal_boot$sd[i])
      f1_pp[i,] <- rnorm(length(f1), mean = f1_boot$mu[i], 
                         sd = f1_boot$sd[i])
      f2_pp[i,] <- rnorm(length(f2), mean = f2_boot$mu[i], 
                         sd = f2_boot$sd[i])
    }else{
      if(lik=="gamma"){
        lon_pp[i,] <- rgamma(length(lon), shape = (lon_boot$mu[i]^2)/(lon_boot$sd[i]^2), 
                             rate = lon_boot$mu[i]/lon_boot$sd[i]^2)
        cal_pp[i,] <- rgamma(length(cal), shape = (cal_boot$mu[i]^2)/(cal_boot$sd[i]^2), 
                             rate = cal_boot$mu[i]/cal_boot$sd[i]^2)
        f1_pp[i,] <- rgamma(length(f1), shape = (f1_boot$mu[i]^2)/(f1_boot$sd[i]^2), 
                            rate = f1_boot$mu[i]/f1_boot$sd[i]^2)
        f2_pp[i,] <- rgamma(length(f2), shape = (f2_boot$mu[i]^2)/(f2_boot$sd[i]^2), 
                            rate = f2_boot$mu[i]/f2_boot$sd[i]^2)
      }else{
        if(lik=="log-normal"){
        lon_pp[i,] <- rlnorm(length(lon), meanlog = log(lon_boot$mu[i]) - .5*log((lon_boot$mu[i]^2 + lon_boot$sd[i]^2)/lon_boot$mu[i]^2),
                             sdlog = sqrt(log((lon_boot$mu[i]^2 + lon_boot$sd[i]^2)/lon_boot$mu[i]^2)))
        cal_pp[i,] <- rlnorm(length(cal), meanlog = log(cal_boot$mu[i]) - .5*log((cal_boot$mu[i]^2 + cal_boot$sd[i]^2)/cal_boot$mu[i]^2),
                             sdlog = sqrt(log((cal_boot$mu[i]^2 + cal_boot$sd[i]^2)/cal_boot$mu[i]^2)))
        f1_pp[i,] <- rlnorm(length(f1), meanlog = log(f1_boot$mu[i]) - .5*log((f1_boot$mu[i]^2 + f1_boot$sd[i]^2)/f1_boot$mu[i]^2),
                            sdlog = sqrt(log((f1_boot$mu[i]^2 + f1_boot$sd[i]^2)/f1_boot$mu[i]^2)))
        f2_pp[i,] <- rlnorm(length(f2), meanlog = log(f2_boot$mu[i]) - .5*log((f2_boot$mu[i]^2 + f2_boot$sd[i]^2)/f2_boot$mu[i]^2),
                            sdlog = sqrt(log((f2_boot$mu[i]^2 + f2_boot$sd[i]^2)/f2_boot$mu[i]^2)))
        }
      } 
    }
    
    delta_pp[i] <- mean(f2_pp[i,]) - (.25*(mean(lon_pp[i,]) + mean(cal_pp[i,])) + .5*mean(f1_pp[i,]))
  }
  
  return(list(lon_pp = lon_pp, cal_pp = cal_pp, f1_pp = f1_pp, f2_pp = f2_pp,
              delta_pp = delta_pp,
              lon = lon,
              cal = cal,
              f1 = f1,
              f2 = f2,
              boots = data.frame(lon = lon_boot$mu,
                                 cal = cal_boot$mu,
                                 f1 = f1_boot$mu,
                                 f2 = f2_boot$mu,
                                 e_f2 = e_f2_boot)))
}


## pp predictive check plots

post_pred_check <- function(df, trait, lik){
  d <- post_pred(df, trait, lik)
  
  lon <- d$lon_pp
  cal <- d$cal_pp
  f1 <- d$f1_pp
  f2 <- d$f2_pp
  
  mu <- sigma <- maximum <- minimum <- matrix(ncol = 4, nrow = 1e4)
  line <- list(lon, cal, f1, f2)
  for(i in 1:4){
    mu[,i] <- apply(line[[i]], 1, mean)
    sigma[,i] <- apply(line[[i]], 1, sd)
    maximum[,i] <- apply(line[[i]], 1, max)
    minimum[,i] <- apply(line[[i]], 1, min)
  }
  
  mu <- data.frame(mu)
  colnames(mu) <- c("lon", "cal", "f1", "f2")
  
  mu_p <- mu %>% 
    pivot_longer(1:4, names_to = "line", values_to = "pp") %>% 
    mutate(ppp = case_when(line=="cal"&pp>=mean(d$cal) ~ 1,
                            line=="lon"&pp>=mean(d$lon) ~ 1,
                            line=="f1"&pp>=mean(d$f1) ~ 1,
                            line=="f2"&pp>=mean(d$f2) ~ 1,
                            T ~ 0)) %>% 
    group_by(line) %>% 
    summarise(upr = quantile(pp, .975),
              lwr = quantile(pp, .025),
              upr.5 = quantile(pp, .75),
              lwr.5 = quantile(pp, .25),
              ppp = as.character(round(sum(ppp)/1e4,2))) %>% 
    mutate(dat = c(mean(d$cal), mean(d$f1), mean(d$f2), mean(d$lon))) %>% 
    ggplot(aes(x = line, y = dat)) +
    geom_point(size = 3) +
    geom_errorbar(aes(x = line, ymax = upr, ymin = lwr), width = 0) +
    geom_errorbar(aes(x = line, ymax = upr.5, ymin = lwr.5), linewidth = 1, width = 0) +
    geom_text(aes(x = line, y = max(upr)+.02, label = ppp)) + 
    theme_classic() +
    labs(x = "Line", 
         title = "Post. Pred. of the Mean", 
         y = "")
  
  sigma <- sigma <- data.frame(sigma)
  colnames(sigma) <- c("lon", "cal", "f1", "f2")
  
  sd_p <- sigma %>% 
    pivot_longer(1:4, names_to = "line", values_to = "pp") %>%  
    mutate(ppp = case_when(line=="cal"&pp>=sd(d$cal) ~ 1,
                           line=="lon"&pp>=sd(d$lon) ~ 1,
                           line=="f1"&pp>=sd(d$f1) ~ 1,
                           line=="f2"&pp>=sd(d$f2) ~ 1,
                           T ~ 0)) %>% 
    group_by(line) %>% 
    summarise(upr = quantile(pp, .975),
              lwr = quantile(pp, .025),
              upr.5 = quantile(pp, .75),
              lwr.5 = quantile(pp, .25),
              ppp = as.character(round(sum(ppp)/1e4,2))) %>% 
    mutate(dat = c(sd(d$cal), sd(d$f1), sd(d$f2), sd(d$lon))) %>% 
    ggplot(aes(x = line, y = dat)) +
    geom_point(size = 3) +
    geom_errorbar(aes(x = line, ymax = upr, ymin = lwr), width = 0) +
    geom_errorbar(aes(x = line, ymax = upr.5, ymin = lwr.5), linewidth = 1, width = 0) +
    geom_text(aes(x = line, y = max(upr)+.02, label = ppp)) + 
    theme_classic() +
    labs(x = "Line", 
         title = "Post. Pred. of the Std.Dev.",
         y = "")
  
  maximum <- data.frame(maximum)
  colnames(maximum) <- c("lon", "cal", "f1", "f2")
  
  max_p <- maximum %>% 
    pivot_longer(1:4, names_to = "line", values_to = "pp") %>%   
    mutate(ppp = case_when(line=="cal"&pp>=max(d$cal) ~ 1,
                           line=="lon"&pp>=max(d$lon) ~ 1,
                           line=="f1"&pp>=max(d$f1) ~ 1,
                           line=="f2"&pp>=max(d$f2) ~ 1,
                           T ~ 0)) %>% 
    group_by(line) %>% 
    summarise(upr = quantile(pp, .975),
              lwr = quantile(pp, .025),
              upr.5 = quantile(pp, .75),
              lwr.5 = quantile(pp, .25),
              ppp = as.character(round(sum(ppp)/1e4,2))) %>%  
    mutate(dat = c(max(d$cal), max(d$f1), max(d$f2), max(d$lon))) %>% 
    ggplot(aes(x = line, y = dat)) +
    geom_point(size = 3) +
    geom_errorbar(aes(x = line, ymax = upr, ymin = lwr), width = 0) +
    geom_errorbar(aes(x = line, ymax = upr.5, ymin = lwr.5), linewidth = 1, width = 0) +
    geom_text(aes(x = line, y = max(upr)+.02, label = ppp)) + 
    theme_classic() +
    labs(x = "Line", 
         title = "Post. Pred. of the Maximum",
         y = "")
  
  minimum <- data.frame(minimum)
  colnames(minimum) <- c("lon", "cal", "f1", "f2")
  
  min_p <- minimum %>% 
    pivot_longer(1:4, names_to = "line", values_to = "pp") %>%   
    mutate(ppp = case_when(line=="cal"&pp>=min(d$cal) ~ 1,
                           line=="lon"&pp>=min(d$lon) ~ 1,
                           line=="f1"&pp>=min(d$f1) ~ 1,
                           line=="f2"&pp>=min(d$f2) ~ 1,
                           T ~ 0)) %>% 
    group_by(line) %>% 
    summarise(upr = quantile(pp, .975),
              lwr = quantile(pp, .025),
              upr.5 = quantile(pp, .75),
              lwr.5 = quantile(pp, .25),
              ppp = as.character(round(sum(ppp)/1e4,2))) %>%
    mutate(dat = c(min(d$cal), min(d$f1), min(d$f2), min(d$lon))) %>% 
    ggplot(aes(x = line, y = dat)) +
    geom_point(size = 3) +
    geom_errorbar(aes(x = line, ymax = upr, ymin = lwr), width = 0) +
    geom_errorbar(aes(x = line, ymax = upr.5, ymin = lwr.5), linewidth = 1, width = 0) +
    geom_text(aes(x = line, y = max(upr)+.02, label = ppp)) + 
    theme_classic() +
    labs(x = "Line", 
         title = "Post. Pred. of the Minimum",
         y = "")
  
  cowplot::plot_grid(mu_p, sd_p, max_p, min_p, ncol = 2)
}

