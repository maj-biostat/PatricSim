

##
# In this design we have x parallel groups


library(configr)
library(optparse)
library(randomizr)
library(doParallel)
library(foreach)
library(dplyr)
library(ggplot2)
library(rstan)
library(drc)
options(mc.cores = 1)
# options(mc.cores = parallel::detectCores()-1)
rstan_options(auto_write = TRUE)
Sys.setenv(LOCAL_CPPFLAGS = '-march=native')

source("setup.R")

tg_env <- new.env()
tg_env$model_code <- rstan::stan_model(file = "logistic_04.stan", auto_write = TRUE)

# tg_env$model_code <- rstan::stan_model(file = "logistic_04b.stan", auto_write = TRUE)
print_tg_env <- function(){
  list(trtgrps = tg_env$trtgrps, 
       lwr = tg_env$lwr,
       upr = tg_env$upr,
       slope = tg_env$slope
       
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




drc_loglogistic <- function(dose = 1, slope = 0, lwr = 0.4, upr = 0.6, ed50 = 5){
  numer = upr - lwr
  
  # when dose is 0, log(dose is -Inf) and 0 times -Inf is NaN so 
  # need to replace with sensible value
  # if slope is zero, denom will be 2
  denom = 1 + exp(slope * (log(dose) - log(ed50)))
  if(slope == 0 & any(is.nan(denom))){
    denom[is.nan(denom)] <- 2
  }
  
  # if slope is zero, response will be mid way between lwr and upr
  res = lwr + (numer / denom)
  as.numeric(res)
}



# library(drc)
# dose = 0:10; slope = 0; lwr = 0.4;  upr = 0.6
# drc_loglogistic(dose, slope, lwr, upr)
# resp <- drc_loglogistic(dose, slope, lwr, upr)
# m <- drm(resp ~ dose, fct = LL.4())
# summary(m)
# plot(dose, drc_loglogistic(dose, slope, lwr, upr), ylim = c(0, 1))

scenario <- function(idx = 1, dose = c(0, 2, 3, 5, 7), 
                     slope = 1, lwr = 0.2, upr = 0.8, ed50 = 5){

  tg_env$slope <- slope
  tg_env$lwr <- lwr
  tg_env$upr <- upr
  
  tg_env$trtgrps <- tibble(
    prob_best = rep(1/length(dose), length(dose)),
    true_mu = drc_loglogistic(dose, slope, lwr, upr),
    prop_rescue = rep(0, length(dose)),
    dose = dose,
    dose_idx = 1:length(dose),
    dose_lab = factor(paste0("D", dose))
  )

  
}

scenario2 <- function(idx = 1, dose = c(0, 2, 3, 5, 7), 
                      location = 4, scale = 2){
  
  tg_env$slope <- slope
  tg_env$lwr <- lwr
  tg_env$upr <- upr
  
  tg_env$trtgrps <- tibble(
    prob_best = rep(1/length(dose), length(dose)),
    true_mu = plogis(dose, location, scale),
    prop_rescue = rep(0, length(dose)),
    dose = dose,
    dose_idx = 1:length(dose),
    dose_lab = factor(paste0("D", dose))
  )
  
  
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



fit_stan <- function(){
  
  # idx = 1; dose = c(1, 2, 3, 5, 7, 10); slope = -3; lwr = 0.4; upr = 0.8; ed50 = 5; scenario(idx, dose, slope, lwr, upr)

  
  print_tg_env()
  
  plot_tg_env_drc()
  
  # redo <- function(s){ 
  #   idx = 1; dose = c(1, 2, 3, 5, 7, 10); slope = s; lwr = 0.4; upr = 0.8; ed50 = 5; scenario(idx, dose, slope, lwr, upr)
  #   
  #   message(print_tg_env())
  #   
  #   print(plot_tg_env_drc())
  # }
  # 
  # redo(-4)
  # 
  # tg_env$df <- generate_trial_data(n_per_arm = 100)
  # gmodels::CrossTable(tg_env$df$y, tg_env$df$dose)
  # tg_env$model_code <- rstan::stan_model(file = "logistic_04b.stan", auto_write = TRUE)
  
  
  # participant data
  # tg_env$df <- generate_trial_data(n_per_arm = 200)
  tmp <- tg_env$df %>%
    dplyr::group_by(dose) %>%
    dplyr::summarise(y = sum(y),
                     trials = n(),
                     prop = y/trials) %>%
    dplyr::ungroup() 
  
  model_data <- list(D = length(tmp$dose),
                     n = tmp$trials,
                     y = tmp$y,
                     tau0 = 2.5)

  # brms::make_stancode(y|trials(trials) ~ dose, data = tmp, family = binomial())
  # grp_means = c(0.2, 0.3, 0.8, 0.4)
  # dt <- tibble(
  #   id = 1:100,
  #   grp = sample(1:4, 100, replace = T)
  # )
  # dt$y <- rnorm(100, mean = grp_means[dt$grp], sd = 1)
  # 
  # brms::make_standata(y ~ 1 + (1|grp), data = dt)
  # brms::make_stancode(y ~ 1 + (1|grp), data = dt)
  
  # model_data <- list(y = tmp$y / tmp$trials,
  #                    trials = tmp$trials,
  #                    dose = tmp %>% dplyr::mutate(dose = exp(dose)) %>% dplyr::pull(dose),
  #                    N = nrow(tmp))
  
  model_fit <- rstan::sampling(tg_env$model_code, 
                               data = model_data,
                               chains = 1, 
                               iter = 10000,
                               refresh = 10000,
                               control = list(adapt_delta = 0.99),
                               verbose = T)
  
  # library(shinystan)
  # shinystan::launch_shinystan(model_fit)
  # print(model_fit, digits = 3)
  # print_tg_env()
  # plot(model_fit, plotfun = "stan_trace", pars = "eta")
  # pairs(model_fit, pars = c("eta"))
  # pairs(model_fit, pars = "yhat")
  
  # log odds of being better by day 7
  # model_draws <- as.matrix(model_fit, pars = c("eta_star"))
  # lapply(1:5, function(x) var(model_draws[,x]))
  # plot_draws <- function(x){
  #   lines(density(model_draws[,x]), col = x + 2)
  # }
  # plot(density(model_draws[,2]), xlim = c(-1.5, 1.5), ylim = c(0, 6))
  # lapply(3:5, plot_draws)
  # model_draws <- as.matrix(model_fit, pars = c("eta"))
  # head(model_draws)
  # plot(density(model_draws[,1]), xlim = c(0, 2.5), ylim = c(0, 5))
  # lapply(2:5, plot_draws)
  
  plot(tmp$dose, tmp$y/tmp$trials, ylim = c(0,1))
  lines(tmp$dose, apply(plogis(model_draws), 2, mean), col = "red")
  lines(tmp$dose, apply(plogis(model_draws), 2, quantile, 0.1), col = "red", lty = 2)
  lines(tmp$dose, apply(plogis(model_draws), 2, quantile, 0.9), col = "red", lty = 2)

}



some_plots <- function(){
  par(mfrow = c(2, 2))
  print(lapply(1:4, function(x) hist(model_draws[, x], main = colnames(model_draws)[x])))
  par(mfrow = c(1, 1))
  
  
  resp_est <- function(dose){
    
    resp_var <- function(x){
      drc_loglogistic(dose, 
                      model_draws[x, "slope"],
                      model_draws[x, "lwr"],
                      model_draws[x, "upr"])
    }
    
    resp <- unlist(lapply(1:nrow(model_draws), resp_var))
    
    # dosex <- rep(dose, length(resp)) 
    # plot(jitter(dosex, 1.3), resp, xlim = c(-0.5, 0.5))
    
    resp
  }
  
  # model_draws_test <- model_draws[sample(1:nrow(model_draws), size = 100), ]
  dfig <- as.data.frame(t(do.call(rbind, lapply(tg_env$trtgrps$dose, resp_est))))
  names(dfig) <- as.character(tg_env$trtgrps$dose_lab)
  dfig1 <- dfig %>%
    tidyr::gather("dose", "response") %>%
    dplyr::mutate(x = substr(dose, 2, nchar(dose)))
  dfig2 <- dfig %>%
    tidyr::gather("dose", "response") %>%
    dplyr::group_by(dose) %>%
    dplyr::summarise(mu = mean(response)) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(x = substr(dose, 2, nchar(dose)))
  
  ggplot(dfig1, aes(x = x, y = response))+
    geom_point() +
    geom_point(data = dfig2, aes(y = mu), col = "red")
}



precision <- function(){
  
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



simulate_trial <- function(id_trial = 1){
  
  # id_trial = 1
  
  message("")
  message("")
  
  message(paste0("###########################################"))
  message(paste0("  TRIAL ", id_trial, " of ", cfg$nsims, " SCENARIO ", cfg$scenarioid))
  message(paste0("###########################################"))
  
  # scenario(cfg$scenarioid)
  scenario(cfg$scenarioid, dose = c(0, 2, 3, 5,7), 
           slope = 1, lwr = 0.2, upr = 0.8)

  # reset data
  tg_env$df <- generate_trial_data()
  
  # gmodels::CrossTable(tg_env$df$rescue_sought, tg_env$df$dose)
  # gmodels::CrossTable(tg_env$df$y, tg_env$df$dose)
  
  # tg_env$trtgrps$n <- tg_env$df %>%
  #   dplyr::group_by(dose) %>%
  #   dplyr::summarise(n = n()) %>%
  #   dplyr::ungroup() %>%
  #   # if one of the groups has not been randomised any patients
  #   # then explicitly set it to zero rather than having it missing
  #   dplyr::right_join(tibble(dose = tg_env$trtgrps$dose), by = "dose") %>%
  #   dplyr::mutate(n = ifelse(is.na(n), 0, n)) %>%
  #   dplyr::pull(n)
  
  fit_stan()
  
 
  tibble(
    id = id_trial, 
    mean_a_diff = tg_env$mean_a_diff,
    mean_p_diff = tg_env$mean_p_diff,
    mean_diff   = tg_env$mean_diff,
    
    prob_a3_ni_a5 = tg_env$prob_a3_ni_a5,
    prob_p3_ni_p5 = tg_env$prob_p3_ni_p5,
    prob_diff_a_eq_diff_p = tg_env$prob_diff_a_eq_diff_p,
    
    a3_ni_a5 = tg_env$a3_ni_a5,
    p3_ni_p5 = tg_env$p3_ni_p5,
    a_eq_p = tg_env$a_eq_p,
    diff_a_eq_diff_p = tg_env$diff_a_eq_diff_p
  )

}




# main loop
# 

interim_seed <- 34432

# myseed <- 1
cfg <- get_cfg()

starttime <- Sys.time()
# set.seed(myseed)

# prealloc list of length nsims
res <- vector(mode = "list", length = cfg$nsims) 

for(i in 1:cfg$nsims){
  
  res[[i]] <- simulate_trial(i)
  
}

dres <- do.call(rbind, res)

endtime <- Sys.time()
difftime(endtime, starttime, units = "hours")
# beepr::beep()
w <- warnings()

rdsfilename <- file.path("out",
                         paste0("design-02-scenario-", cfg$scenarioid, "-",
                                format(Sys.time(), "%Y-%m-%d-%H-%M-%S"), ".RDS"))

saveRDS(list(scenario = cfg$scenarioid,
             results = dres,
             trtgrps = tg_env$trtgrps,
             warnings = w,
             starttime = starttime,
             endtime = endtime,
             duration = difftime(endtime, starttime, units = "hours")),
        rdsfilename)
assign("last.warning", NULL, envir = baseenv())

  
  
