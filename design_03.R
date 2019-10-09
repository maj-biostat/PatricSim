

##
# In this design we have x parallel groups

# https://www.researchgate.net/post/dose_response_curves_is_better_to_create_different_models_or_just_one
# DRC: make sure you are using the right function, as the drc
# can provide a number of drc_loglogistical fits, ie. Weibull, logistic, drc_loglogistic.

# The next step is to describe your data in the most parsimonious way -  this
# means that you are required to do a model reduction. In your example you have
# shown a five parameter drc_loglogistic function, which will have to be compared
# with a four, a three and perhaps also a two parameter model.

# When reducing your model, various parameters become 0. In the instance of your
# four parameter model, you have the upper limit (=d), lower limit (=c), slope
# (=b), and the inception point (=e, this parameter also is known as the ED50 or
# median point in the curve decay). The three parameter model assumes that your
# lower limit c is equal to 0, and the two parameter model you will only
# quantify the slope and the inception point (this is probably too simple).

# The first step is to compare the various parameter functions using the
# mselect() function (e.g. mselect(m1, list(LL.5(), LL.4(), LL.3(), LL.2()),
# linreg=TRUE, icfct=AIC). This will help you reduce the model to the
# appropriate parameter function, and in the example I have given you, based on
# the AIC, loglikelihood and lack of fit. Also note that the linreg =TRUE means
# that the curves are also compared against linear functions. The function will
# rank your curves based on the best fit. Caution here, as this is user defined,
# so you need to check that you are using a biologically appropriate function.
# You can also use mselect() to compare between different families if you need
# to.

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
tg_env$model_code <- rstan::stan_model(file = "logistic_03.stan", auto_write = TRUE)


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



drc_loglogistic <- function(dose = 1, slope = 4, lwr = 0.2, upr = 0.7){
  numer = upr - lwr
  
  denom = 1 + exp(slope * log(dose))
  if(slope == 0 & any(is.nan(denom))){
    denom[is.nan(denom)] <- 2
  }
  
  res = lwr + (numer / denom)
  as.numeric(res)
}


# dose = 0:10; slope = 0; lwr = 0.2;  upr = 0.8
# drc_loglogistic(dose, slope, lwr, upr)
# plot(dose, drc_loglogistic(dose, slope, lwr, upr), ylim = c(0, 1))

scenario <- function(idx = 1, dose = c(0, 2, 3, 5,7), 
                     slope = 1, lwr = 0.2, upr = 0.8){
  
  
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
  
  # idx = 1; dose = c(0, 2, 3, 5, 7); slope = -1; lwr = 0.7; upr = 0.8; scenario(idx, dose, slope, lwr, upr)

  
  print_tg_env()
  
  plot_tg_env_drc()
  
  # tg_env$df <- generate_trial_data(n_per_arm = 100)
  # gmodels::CrossTable(tg_env$df$y, tg_env$df$dose)
  # tg_env$model_code <- rstan::stan_model(file = "logistic_03.stan", auto_write = TRUE)
  
  
  # participant data
  # tg_env$df <- generate_trial_data(n_per_arm = 200)
  tmp <- tg_env$df %>%
    dplyr::group_by(dose) %>%
    dplyr::summarise(y = sum(y),
                     trials = n(),
                     prop = y/trials) %>%
    dplyr::ungroup() 
  
  model_data <- list(y = tmp$y ,
                     trials = tmp$trials,
                     dose = tmp %>%  dplyr::pull(dose),
                     N = nrow(tmp))
  
  
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
  
  # print(model_fit, digits = 3)
  # print_tg_env()
  # plot(model_fit, plotfun = "stan_trace", pars = c("slope", "lwr", "upr"))
  # pairs(model_fit, pars = c("slope", "lwr", "upr"))
  # pairs(model_fit, pars = "yhat")
  
  # log odds of being better by day 7
  model_draws <- as.matrix(model_fit, pars = c("slope", "lwr", "upr"))
  head(model_draws)
  
  
  
  
  

}



some_plots <- function(){
  par(mfrow = c(2, 2))
  print(lapply(1:3, function(x) hist(model_draws[, x], main = colnames(model_draws)[x])))
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

  
  
