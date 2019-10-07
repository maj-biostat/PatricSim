

##
# In this design we have x parallel groups


library(configr)
library(optparse)
library(randomizr)
library(dplyr)
library(ggplot2)
library(rstan)
library(walker)
#library(drc)
options(mc.cores = 1)
# options(mc.cores = parallel::detectCores()-1)
rstan_options(auto_write = TRUE)
Sys.setenv(LOCAL_CPPFLAGS = '-march=native')

source("setup.R")

tg_env <- new.env()
# tg_env$model_code <- rstan::stan_model(file = "design_04_3.stan", auto_write = TRUE)




precision_test <- function(){
  
  # 
  nsim <- 1000
  
  # assume a uniform prior
  a <- 1
  b <- 1
  
  # true recovery 80%
  p <- 0.8
  n <- 100
  
  m <- array(0, dim = c(nsim, 2))
  
  for(i in 1:nsim){
    evt <- rbinom(1, n, p)
    
    y <- rbeta(1000, a + evt, b + n - evt)
    
    m[i, ] <- quantile(y, probs = c(0.025, 0.975))
  }
  
  
  ci_est <- colMeans(m)
  ci_est
  diff(ci_est)
  
}



print_tg_env <- function(){
  list(trtgrps = tg_env$trtgrps, 
       location = tg_env$location,
       scale = tg_env$scale,
       p_range = tg_env$p_range,
       p_lwr = tg_env$p_lwr
  )
}

plot_tg_env_drc <- function(){
  
  plot(tg_env$trtgrps$dose, tg_env$trtgrps$true_mu, ylim = c(0, 1))
}

estBetaParams <- function(mu, sd) {
  var = sd^2
  
  if(mu < 0.00001 | mu > 0.99999) return(NA)
  if(var < 0.00001 | var > 0.24999) return(NA)
  
  alpha <- ((1 - mu) / var - 1 / mu) * mu ^ 2
  beta <- alpha * (1 / mu - 1)
  return(params = list(alpha = alpha, beta = beta))
}

drc_cumnorm <- function(dose, 
                        location, scale, 
                        p_range, p_lwr, plot = F){
  
  p_rec <- p_range * pnorm(dose, location, scale) + p_lwr
  
  if(plot){
    plot(dose, p_rec, 
         ylim = 0:1, type = "l", 
         xlab = "Duration", ylab = "Pr(recov)")
    points(dose, p_rec)
  }
  
  return(p_rec)
  
}

scenario <- function(idx = 0, dose = 0:10){
  

  tg_env$trtgrps <- tibble(
    prob_best = rep(1/length(dose), length(dose)),
    true_mu = rep(0, length(dose)),
    prop_rescue = rep(0, length(dose)),
    dose = dose,
    dose_idx = 1:length(dose),
    dose_lab = factor(paste0("D", dose))
  )

 
  if(idx == 1){
    tg_env$location <- 5
    tg_env$scale <- 1.2
    tg_env$p_range <- 0
    tg_env$p_lwr <- 0.5
  
  } else if(idx == 2){
    
    tg_env$location <- 3
    tg_env$scale <- 1.2
    tg_env$p_range <- 0.1
    tg_env$p_lwr <- 0.5

  } else if(idx == 3){
    
    tg_env$location <- 5
    tg_env$scale <- 1.2
    tg_env$p_range <- 0.1
    tg_env$p_lwr <- 0.5

  } else if(idx == 4){
    
    tg_env$location <- 7
    tg_env$scale <- 1.2
    tg_env$p_range <- 0.1
    tg_env$p_lwr <- 0.5

  } else if(idx == 5){
    
    tg_env$location <- 3
    tg_env$scale <- 1.2
    tg_env$p_range <- 0.3
    tg_env$p_lwr <- 0.5

  } else if(idx == 6){
    
    tg_env$location <- 5
    tg_env$scale <- 1.2
    tg_env$p_range <- 0.3
    tg_env$p_lwr <- 0.5
    
  } else if(idx == 7){
    
    tg_env$location <- 7
    tg_env$scale <- 1.2
    tg_env$p_range <- 0.3
    tg_env$p_lwr <- 0.5
    
  } 
  
  tg_env$trtgrps$true_mu <- drc_cumnorm(dose = tg_env$trtgrps$dose,
                                        location = tg_env$location, 
                                        scale = tg_env$scale, 
                                        p_range = tg_env$p_range, 
                                        p_lwr = tg_env$p_lwr)
}



generate_trial_data <- function(n_per_arm = 100) {
  
  
  n <- nrow(tg_env$trtgrps)*n_per_arm
  dat <- tibble(
    id = 1:n
  )
  
  dat$grp <- randomizr::complete_ra(N = n, 
                                    num_arms = nrow(tg_env$trtgrps),
                                    conditions = tg_env$trtgrps$dose_lab)
  
  dat$dose_idx <- as.numeric(dat$grp)
  
  dat$dose <- tg_env$trtgrps$dose[dat$dose_idx]
  dat$y = rbinom(n = n, size = 1, prob = tg_env$trtgrps$true_mu[dat$dose_idx])
  dat$prop_rescue <- tg_env$trtgrps$prop_rescue[dat$dose_idx]
  
  idx_failures <- dat %>%
    dplyr::group_by(grp) %>%
    dplyr::sample_frac(size = prop_rescue) %>%
    dplyr::ungroup() %>%
    dplyr::pull(id)
  # split by trt group
  # sample prop_rescue records and set y to zero
  
  dat$rescue_sought <- 0
  dat$rescue_sought[dat$id %in% idx_failures] <- 1
  
  # gmodels::CrossTable(dat$rescue_sought, dat$grp)
  
  dat$y_orig <- dat$y
  dat$y[dat$id %in% idx_failures] <- 0
  dat
}


p_best <- function(mat) {
  as.numeric(prop.table(table(factor(max.col(mat), levels = 1:ncol(mat)))))
}


# Tests bias and coverage.
model_test <- function(){
  
  get_ci <- function(x){
    
    tg_env$df <- generate_trial_data(n_per_arm = 100)
    
    tmp <- tg_env$df %>%
      dplyr::group_by(dose) %>%
      dplyr::summarise(y = sum(y),
                       trials = n(),
                       prop = y/trials) %>%
      dplyr::ungroup() 
    
    model_data <- list(D = length(tmp$dose),
                       duration = tmp$dose,
                       n = tmp$trials,
                       y = tmp$y,
                       tau0 = 2.5)
    
    model_fit <- rstan::sampling(tg_env$model_code, 
                                 data = model_data,
                                 chains = 1, 
                                 iter = cfg$mcmc_iter,
                                 refresh = cfg$mcmc_iter,
                                 control = list(adapt_delta = 0.99),
                                 verbose = T)
    
    model_draws <- as.matrix(model_fit, pars = c("eta"))
    
    mu <- apply(plogis(model_draws), 2, mean)
    
    ci_mu <- rbind(apply(plogis(model_draws), 2, quantile, 0.05),
                   apply(plogis(model_draws), 2, quantile, 0.975))
    
    ci_cov <- ci_mu[1, ] < mu & ci_mu[2, ] > mu
    
    bias <- tg_env$trtgrps$true_mu - mu
    
    rbind(ci_cov, bias)
    
  }
  
  res <- do.call(rbind, lapply(1:100, get_ci))
  
  
  
}

fit_stan_1 <- function(){

  # participant data
  tmp <- tg_env$df %>%
    dplyr::group_by(dose) %>%
    dplyr::summarise(y = sum(y),
                     trials = n(),
                     prop = y/trials) %>%
    dplyr::ungroup() 
  
 #  plot(tmp$dose, tmp$prop, ylim = 0:1)
  
  model_data <- list(D = length(tmp$dose),
                     duration = tmp$dose,
                     n = tmp$trials,
                     y = tmp$y,
                     tau0 = 2.5)

  model_fit <- rstan::sampling(tg_env$model_code, 
                               data = model_data,
                               chains = 1, 
                               iter = cfg$mcmc_iter,
                               refresh = cfg$mcmc_iter,
                               control = list(adapt_delta = 0.99),
                               verbose = T)
  
  # library(shinystan)
  # shinystan::launch_shinystan(model_fit)
  # print(model_fit, digits = 3)
  # print_tg_env()
  # plot(model_fit, plotfun = "stan_trace", pars = "eta")
  # pairs(model_fit, pars = c("eta"))
  
  # log odds of being better by day 7
  prop_recov <- plogis(as.matrix(model_fit, pars = c("eta")))
  
  prob_recov_mu <- colMeans(prop_recov)
  
  # Differences between proportions recovered - 
  # computes T_dmax - T_d with d < dmax
  prop_recov_diffs <- prop_recov[, ncol(prop_recov)] - prop_recov[, 1:(ncol(prop_recov)-1)]
  
  # superiority 
  # NOTE!! ordered as p_11 - p1, p11 - p2, p_11 - p_3 etc
  prob_sup <- colMeans(prop_recov_diffs > 0)
  # decis_sup <- colMeans(prop_recov_diffs > 0) > tg_env$p_sup_deicsion_thresh
  
  # NOTE!! ordered as p_1 ni p11, p_2 ni p11, p_3 ni p11, 
  prob_ni <- colMeans(prop_recov_diffs > -cfg$p_ni_thresh & 
                        prop_recov_diffs < cfg$p_ni_thresh)
  
  # decis_ni <- colMeans(prop_recov_diffs > -tg_env$p_ni_thresh & 
  #            prop_recov_diffs < tg_env$p_ni_thresh) > tg_env$p_ni_deicsion_thresh
  
  dres <- rbind(prob_sup, prob_ni)
  dres
}


fit_stan_2 <- function(){
  
  # participant data
  tmp <- tg_env$df %>%
    dplyr::group_by(dose) %>%
    dplyr::summarise(y = sum(y),
                     trials = n(),
                     prop = y/trials) %>%
    dplyr::ungroup() 
  
  #  plot(tmp$dose, tmp$prop, ylim = 0:1)
  
  model_data <- list(D = length(tmp$dose),
                     duration = tmp$dose,
                     n = tmp$trials,
                     y = tmp$y)
  
  # tg_env$model_code <- rstan::stan_model(file = "design_04_2.stan", auto_write = TRUE)
  model_fit <- rstan::sampling(tg_env$model_code, 
                               data = model_data,
                               chains = 1, 
                               iter = cfg$mcmc_iter,
                               refresh = cfg$mcmc_iter,
                               control = list(adapt_delta = 0.99),
                               verbose = T)
  
  # library(shinystan)
  # shinystan::launch_shinystan(model_fit)
  print(model_fit, digits = 3)
  # print_tg_env()
  # plot(model_fit, plotfun = "stan_trace", pars = "eta")
  # pairs(model_fit, pars = c("theta_1"))
  
  # log odds of being better by day 7
  prop_recov <- plogis(as.matrix(model_fit, pars = c("theta")))

  prob_recov_mu <- colMeans(prop_recov)
  prob_recov_mu_lwr <- apply(prop_recov, 2, quantile, 0.1)
  prob_recov_mu_upr <- apply(prop_recov, 2, quantile, 0.9)
  
  bias_mu <- prob_recov_mu - tg_env$trtgrps$true_mu
  
  # Differences between proportions recovered - 
  # computes T_dmax - T_d with d < dmax
  prop_recov_diffs <- prop_recov[, ncol(prop_recov)] - prop_recov[, 1:(ncol(prop_recov)-1)]
  
  # superiority 
  # NOTE!! ordered as p_11 - p1, p11 - p2, p_11 - p_3 etc
  prob_sup <- colMeans(prop_recov_diffs > 0)
  prob_sup
  # decis_sup <- colMeans(prop_recov_diffs > 0) > tg_env$p_sup_deicsion_thresh
  
  # NOTE!! ordered as p_1 ni p11, p_2 ni p11, p_3 ni p11, 
  prob_ni <- colMeans(prop_recov_diffs > -cfg$p_ni_thresh & 
                        prop_recov_diffs < cfg$p_ni_thresh)
  prob_ni
  # decis_ni <- colMeans(prop_recov_diffs > -tg_env$p_ni_thresh & 
  #            prop_recov_diffs < tg_env$p_ni_thresh) > tg_env$p_ni_deicsion_thresh
  
  dres <- rbind(prob_sup, prob_ni)
  dres
  # plot(tg_env$trtgrps$dose, tmp$prop, ylim = c(0,1))
  # lines(tg_env$trtgrps$dose, tg_env$trtgrps$true_mu)
  # lines(tg_env$trtgrps$dose, apply(prop_recov, 2, mean), col = "red")
  # lines(tg_env$trtgrps$dose, apply(prop_recov, 2, quantile, 0.1), col = "red", lty = 2)
  # lines(tg_env$trtgrps$dose, apply(prop_recov, 2, quantile, 0.9), col = "red", lty = 2)
  
  
  return(list(dres = dres,
              prob_recov_mu = prob_recov_mu, 
              prob_recov_mu_lwr = prob_recov_mu_lwr, 
              prob_recov_mu_upr = prob_recov_mu_upr, 
              bias_mu = bias_mu))
}



fit_walker_1 <- function(){
  
  # participant data
  tmp <- tg_env$df %>%
    dplyr::group_by(dose) %>%
    dplyr::summarise(y = sum(y),
                     trials = n(),
                     prop = y/trials) %>%
    dplyr::ungroup() 
  
  #  plot(tmp$dose, tmp$prop, ylim = 0:1)
  
  model_data <- list(D = length(tmp$dose),
                     duration = tmp$dose,
                     n = tmp$trials,
                     y = tmp$y)
  
  # tg_env$model_code <- rstan::stan_model(file = "design_04_2.stan", auto_write = TRUE)
  model_fit <- walker_glm(y ~ -1 + 
                            rw1(~ dose, beta_prior = c(0, 10), sigma_prior = c(0, 10)), 
                          data = tmp, u = tmp$trials, distribution = "binomial",
                          refresh = 0, chains = 1)
  
  # plot_fit(model_fit)
  # print(model_fit$stanfit) 
  
  # library(shinystan)
  # shinystan::launch_shinystan(model_fit$stanfit)
  # print(model_fit, digits = 3)
  # print_tg_env()
  # plot(model_fit, plotfun = "stan_trace", pars = "eta")
  # pairs(model_fit, pars = c("theta_1"))
  
  # log odds of being better by day 7
  prop_recov <- as.matrix(model_fit$stanfit, pars = c("y_fit"))

  prob_recov_mu <- colMeans(prop_recov)
  prob_recov_mu_lwr <- apply(prop_recov, 2, quantile, 0.1)
  prob_recov_mu_upr <- apply(prop_recov, 2, quantile, 0.9)
  
  bias_mu <- prob_recov_mu - tg_env$trtgrps$true_mu
  
  # Differences between proportions recovered - 
  # computes T_dmax - T_d with d < dmax
  prop_recov_diffs <- prop_recov[, ncol(prop_recov)] - prop_recov[, 1:(ncol(prop_recov)-1)]
  
  # superiority 
  # NOTE!! ordered as p_11 - p1, p11 - p2, p_11 - p_3 etc
  prob_sup <- colMeans(prop_recov_diffs > 0)
  
  # decis_sup <- colMeans(prop_recov_diffs > 0) > tg_env$p_sup_deicsion_thresh
  
  # NOTE!! ordered as p_1 ni p11, p_2 ni p11, p_3 ni p11, 
  prob_ni <- colMeans(prop_recov_diffs > -cfg$p_ni_thresh & 
                        prop_recov_diffs < cfg$p_ni_thresh)
  
  # decis_ni <- colMeans(prop_recov_diffs > -tg_env$p_ni_thresh & 
  #            prop_recov_diffs < tg_env$p_ni_thresh) > tg_env$p_ni_deicsion_thresh
  
  dres <- rbind(prob_sup, prob_ni)
  
  
  # plot(tg_env$trtgrps$dose, tmp$prop, ylim = c(0,1))
  # lines(tg_env$trtgrps$dose, tg_env$trtgrps$true_mu)
  # lines(tg_env$trtgrps$dose, apply(prop_recov, 2, mean), col = "red")
  # lines(tg_env$trtgrps$dose, apply(prop_recov, 2, quantile, 0.1), col = "red", lty = 2)
  # lines(tg_env$trtgrps$dose, apply(prop_recov, 2, quantile, 0.9), col = "red", lty = 2)
  
  
  return(list(dres = dres,
              prob_recov_mu = prob_recov_mu, 
              prob_recov_mu_lwr = prob_recov_mu_lwr, 
              prob_recov_mu_upr = prob_recov_mu_upr, 
              bias_mu = bias_mu))
}


simulate_trial <- function(id_trial = 1){
  
  # id_trial = 1
  
  message("")
  message("")
  
  message(paste0("###########################################"))
  message(paste0("  TRIAL ", id_trial, " of ", cfg$nsims, " SCENARIO ", cfg$scenarioid))
  message(paste0("###########################################"))
  
  #plot_tg_env_drc()

  # reset data - DEFAULTS TO 100 PER ARM
  tg_env$df <- generate_trial_data(cfg$n_per_trt)

  l <- fit_walker_1()

  l

}




# main loop
# 

# interim_seed <- 34432

# myseed <- 1
cfg <- get_cfg()

scenario(cfg$scenarioid, 
         dose = cfg$dose)

cfg
print_tg_env()

starttime <- Sys.time()
# set.seed(myseed)

# prealloc list of length nsims
res <- vector(mode = "list", length = cfg$nsims) 

for(i in 1:cfg$nsims){
  
  res[[i]] <- simulate_trial(i)
  
}


unpack1 <- function(x){
  res[[x]]$dres
}

dtmp <- do.call(rbind, lapply(1:length(res), unpack1))
dsup <- dtmp[seq(1, nrow(dtmp), by = 2),]
dni <- dtmp[seq(2, nrow(dtmp), by = 2),]

unpack2 <- function(x){
  res[[x]]$prob_recov_mu
}

dprob <- do.call(rbind, lapply(1:length(res), unpack2))

unpack3 <- function(x){
  res[[x]]$prob_recov_mu_lwr
}

dprob_lwr <- do.call(rbind, lapply(1:length(res), unpack3))

unpack4 <- function(x){
  res[[x]]$prob_recov_mu_upr
}

dprob_upr <- do.call(rbind, lapply(1:length(res), unpack4))

unpack5 <- function(x){
  res[[x]]$bias_mu
}

dbias <- do.call(rbind, lapply(1:length(res), unpack5))


endtime <- Sys.time()
difftime(endtime, starttime, units = "hours")
# beepr::beep()
w <- warnings()

rdsfilename <- file.path("out",
                         paste0("design-04-scenario-", cfg$scenarioid, "-",
                                format(Sys.time(), "%Y-%m-%d-%H-%M-%S"), ".RDS"))

saveRDS(list(scenario = cfg$scenarioid,
             dsup = dsup,
             dni = dni,
             dprob = dprob,
             dprob_lwr = dprob_lwr,
             dprob_upr = dprob_upr,
             dbias = dbias,
             trtgrps = list(trtgrps = tg_env$trtgrps,
                            location = tg_env$location,
                            scale = tg_env$scale,
                            p_range = tg_env$p_range,
                            p_lwr = tg_env$p_lwr),
             nsim = cfg$nsims,
             n_per_trt = cfg$n_per_trt,
             p_ni_thresh = cfg$p_ni_thresh,
             p_ni_deicsion_thresh = cfg$p_ni_deicsion_thresh,
             p_sup_deicsion_thresh = cfg$p_sup_deicsion_thresh,
             warnings = w,
             starttime = starttime,
             endtime = endtime,
             duration = difftime(endtime, starttime, units = "hours")),
        rdsfilename)
assign("last.warning", NULL, envir = baseenv())

  
  
