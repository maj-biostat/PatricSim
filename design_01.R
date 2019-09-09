

##
# In this design we have 4 groups representing 7, 5, 3 and zero days
# of treatment. The analysis is based on the two bottle design. 
# 7 days treatment is the reference level that is used in comparisons.
# Randomisation is constant to each arm until a minmimum number of 
# participants is reached in each arm and then randomisation is by RAR.
# The response is binary, representing whether a participant recovered
# by day 7 or not. 
# A simple logistic regression model is used for the analysis. No intercept
# is specified and so the parameter estimates give the group means of 
# log-odds of recovery.
# We compute the probability that each group is the best. We also run a
# non-inferiority test by first computing the proportion recovered
# in each group and then comparing each group to the reference arm of 7 days.
# From this we get the probabilities that each arm is non-inferior to the
# reference arm.
# At each interim we do NI tests. We check whether the proportion recovered
# in each arm is within 0.1 of the reference arm. If the probability of NI 
# is > 0.85 then the arm is considered NI. If we conclude that all 
# arms are NI to the reference arm, then we halt the trial.
# A number of scenarios is defined covering no differences in the true 
# proportion recovered, a superior reference arm, non-inferiority in one of
# the lower dose durations etc.

library(configr)
library(optparse)
library(randomizr)
library(doParallel)
library(foreach)
library(dplyr)
library(ggplot2)
library(rstan)
options(mc.cores = 1)
# options(mc.cores = parallel::detectCores()-1)
rstan_options(auto_write = TRUE)
Sys.setenv(LOCAL_CPPFLAGS = '-march=native')

source("setup.R")

# define treatment groups, initial allocation probability and true means
grp <- 1:4

tg_env <- new.env()

# min samp size for rar by subgroup
rar_min_cell <- 50
non_inf_delta <- 0.1

# tg_env$trtgrps
# tg_env$rar_active
# tmp2 <- tg_env$trtgrps


tg_env$model_code <- rstan::stan_model(file = "logistic.stan", 
                                       auto_write = TRUE)







scenario <- function(idx = 1){
  
  tg_env$trtgrps <- tibble(
    prob_best = rep(1/length(grp), length(grp)),
    rand_prob = rep(1/length(grp), length(grp)),
    true_mean = rep(0.8, length(grp)),
    n = 0,
    grp = as.factor(paste0("T", 1:length(grp))),
    grplab = as.character(paste0(c(7, 5, 3, 0), "-day")),
    # logodds and var used for pbest and rar
    est_logodds = rep(0, length(grp)),
    est_var = rep(0, length(grp)),
    est_prop = rep(0, length(grp)),
    prob_ni = rep(0, length(grp)),
    trt_ni = rep(0, length(grp))
  )
  
  rownames(tg_env$trtgrps) <- NULL
  tg_env$rar_active <- 0
  # we have to have this much certainty that the
  # alt treatments are non-inferior using a non_inf_delta
  tg_env$ni_thresh <- 0.85
  
  
  if(idx == 1){
    # all same - high prob of recovery -
    # 0.9 will give you estimation issues, e.g. all values in grp set to 1
    tg_env$trtgrps$true_mean = rep(0.8, length(grp))
  } else if(idx == 2) {
    # all same - med prob of recovery
    tg_env$trtgrps$true_mean = rep(0.7, length(grp))
  } else if(idx == 3) {
    # all same - low prob of recovery
    tg_env$trtgrps$true_mean = rep(0.6, length(grp))
  } else if(idx == 4) {
    # 7 day ab is better than everything
    tg_env$trtgrps$true_mean = c(0.8, 0.5, 0.5, 0.5)
  } else if(idx == 5) {
    # 7 day ab is no different from 5 day, others are worse
    tg_env$trtgrps$true_mean = c(0.7, 0.7, 0.5, 0.5)
  } else if(idx == 6) {
    # any ab is better than none
    tg_env$trtgrps$true_mean = c(0.7, 0.7, 0.7, 0.5)
  } else if(idx == 7) {
    # any ab is better than none gradient
    tg_env$trtgrps$true_mean = c(0.9, 0.8, 0.7, 0.5)
  }
}



generate_trial_data <- function(idx_start = 1, idx_end = 100) {
  
  
  n <- idx_end - idx_start + 1
  dat <- tibble(
    id = idx_start:idx_end
  )
  
  if(!tg_env$rar_active){
    dat$grp <- randomizr::complete_ra(N = n, num_arms = length(grp))
  } else {
    
    dat$grp <- randomizr::complete_ra(N = n,
                                      num_arms = length(grp),
                                      prob_each = tg_env$trtgrps$rand_prob)
    
    message(paste0("Groups produced by randomizr ", length(unique(dat$grp))))
    
    # if(length(unique(dat$grp)) < max(grp)){
    #   rngdat <- T
    #   its <- 1
    #   while(rngdat){
    #     dat$grp <- randomizr::complete_ra(N = n,
    #                                       num_arms = length(grp),
    #                                       prob_each = tg_env$trtgrps$rand_prob)
    #     if(length(unique(dat$grp)) == max(grp)){
    #       rngdat <- F
    #     }
    #     if(its > 10){
    #       message("dbg")
    #     }
    #     its <- its + 1
    #   }
    # }
    
    
    
  }
  
  # prevents the case where all pts in a given group are all successes
  rngdat <- T
  its <- 1
  grpidx <- as.numeric(dat$grp)
  while(rngdat){
    message(paste0("Data generation attempt ", its))
    dat$y = rbinom(n = n, size = 1, prob = tg_env$trtgrps$true_mean[grpidx])
    
    if(nrow(tg_env$df)>0){
      tmp <- rbind(tg_env$df, dat)
    } else{
      tmp <- dat
    }
    
    cell_n <- tmp %>%
      dplyr::group_by(y, grp) %>%
      dplyr::summarise(n = n()) %>%
      dplyr::ungroup() %>%
      dplyr::pull(n)
    
    message(paste0(" ", cell_n))
    
    if(length(cell_n) == length(grp)*2 & all(cell_n > 0)){
      rngdat <- F
    }
    
    its <- its + 1

    
    if(its > 10){

      saveRDS(list(dat, cell_n, idx_start, idx_end, grpidx,
                   tg_env$trtgrps$true_mean),
              file = "panic.RDS")
      stop(paste0("Problem generating data."))
      
      
    }
  }
  
  # aggregate(dat$y, list(dat$grp), mean)
  # plot(table(dat$grp), ylim = c(0, 60))
  # plot(aggregate(dat$y, list(dat$grp), mean))
  
  dat
}


p_best <- function(mat) {
  as.numeric(prop.table(table(factor(max.col(mat), levels = 1:ncol(mat)))))
}

# inv_logit <- function(x){
#   exp(x)/(1+exp(x))
# }
#
# odd <- function(x){
#   x/(1-x)
# }



fit_stan <- function(){
  
  # participant data
  Xmat <- model.matrix(y ~ 0 + grp, data = tg_env$df)
  
  # lm1 <- glm(y ~ 0 + grp, data = df, family = "binomial")
  # summary(lm1)
  
  ## Zmat for estimating mean change in each group
  Zmat <- as.data.frame(model.matrix(~ 0 + grp, data = tg_env$trtgrps))
  
  model_data <- list(y = tg_env$df$y,
                     X = Xmat,
                     N = length(tg_env$df$y), K = 4, prior_only = 0)
  
  model_fit <- rstan::sampling(tg_env$model_code, 
                               data = model_data,
                               chains = 1, 
                               iter = 2000,
                               refresh = 2000,
                               seed = interim_seed,
                               control = list(adapt_delta = 0.999),
                               verbose = F)
  
  # print(model_fit)
  
  # log odds of being better by day 7
  model_draws <- as.matrix(model_fit, pars = c("b"))
  # next bit of matrix mult redundant in that model_draws is what
  # we want but kept it in the event that new predictors introduced.
  # log odds of being better by day 7
  # model_mu <- model_draws %*% t(Zmat)
  model_mu <- model_draws
  # colMeans(model_mu)
  
  model_props <- apply(model_mu, 2, function(x) exp(x)/(1+exp(x)))
  colnames(model_props) <- paste0("b", 1:4)
  tg_env$trtgrps$est_prop <- colMeans(model_props)
  
  # comparisons
  model_prop_diffs <- do.call(cbind, lapply(1:4, function(x)
    model_props[,1] - model_props[,x]))
  
  tg_env$trtgrps$est_logodds <- apply(model_mu, 2, mean)
  tg_env$trtgrps$est_var <- diag(var(model_mu))
  tg_env$trtgrps$prob_best <- p_best(model_mu)
  tg_env$trtgrps$prob_ni <- colMeans(model_prop_diffs < non_inf_delta)

  if(tg_env$rar_active){
    
    message("Original alloc prob = ", sprintf("%.3f ", tg_env$trtgrps$rand_prob))
    
    tg_env$trtgrps$rand_prob <- sqrt((tg_env$trtgrps$prob_best *
                                        tg_env$trtgrps$est_var)/(tg_env$trtgrps$n + 1))
    tg_env$trtgrps$rand_prob <- tg_env$trtgrps$rand_prob / sum(tg_env$trtgrps$rand_prob)
    
    # if(any(tg_env$trtgrps$rand_prob == 0)){
    #   message("Fixing arm with zero prob ")
    #   tg_env$trtgrps$rand_prob[tg_env$trtgrps$rand_prob == 0] <- 0.01  
    #   tg_env$trtgrps$rand_prob <- tg_env$trtgrps$rand_prob / sum(tg_env$trtgrps$rand_prob)
    # }
    
    message("Updated alloc prob  = ", sprintf("%.3f ", tg_env$trtgrps$rand_prob))
    
  } else {
    message("Rar not active n = ", paste0(tg_env$trtgrps$n, sep = " "))
    
  }
 
}



simulate_trial <- function(id_trial){
  
  # scenarioid <- 1
  
  message(paste0("#################################"))
  message(paste0("TRIAL ", id_trial, " SCENARIO ", cfg$scenarioid))
  message(paste0("#################################"))
  
  scenario(cfg$scenarioid)
  
  looks <- seq(100,500,by=100)
  
  # reset data
  tg_env$df <- tibble()
  
  for(idx_intrm in 1:length(looks)){
    
    # idx_intrm=1
    # idx_intrm=idx_intrm+1
    # idx_intrm=2
    # idx_intrm=3
    message(paste0("Interim ", idx_intrm))
    tg_env$df <- rbind(tg_env$df,
                generate_trial_data(idx_start = nrow(tg_env$df)+1,
                                    idx_end = looks[idx_intrm]))

    
    # rar disabled until each cell has rar_min_cell (10)
    tg_env$rar_active <- all(as.numeric(table(tg_env$df$grp)) > rar_min_cell)
    
    # cumulative total of participants assigned to each group
    tg_env$trtgrps$n <- tg_env$df %>%
      dplyr::group_by(grp) %>%
      dplyr::summarise(n = n()) %>%
      dplyr::ungroup() %>%
      # if one of the groups has not been randomised any patients
      # then explicitly set it to zero rather than having it missing
      dplyr::right_join(tibble(grp = tg_env$trtgrps$grp), by = "grp") %>%
      dplyr::mutate(n = ifelse(is.na(n), 0, n)) %>%
      dplyr::pull(n)
    
    # table(tg_env$df$subX)
    
    fit_stan()
    
    trts_ni <- c(tg_env$trtgrps$prob_ni > tg_env$ni_thresh)
    if(all(trts_ni)){
      tg_env$trtgrps$trt_ni <- trts_ni
      tg_env$trtgrps$trt_ni[1] <- NA
      message("All treatments NI: ", paste0(trts_ni, sep= " "))
      message("Exiting current trial simulation.")
      break
    }
    
  }
  
  
  cbind(id = id_trial, tg_env$trtgrps)
  
  
}

figs <- function(){
  dfig <- tg_env$df %>%
    dplyr::mutate(grplab = tg_env$trtgrps$grplab[as.numeric(gsub("T", "", grp))]) %>%
    dplyr::group_by(grp, grplab) %>%
    dplyr::summarise(n = n(),
                     yhat = mean(y)) %>%
    dplyr::ungroup()
  
  ggplot(dfig, aes(x = grp, y = n)) +
    geom_point() +
    scale_y_continuous("n", limits = c(0, max(tg_env$trtgrps$n)))+
    ggtitle("Samples Size")
  
  ggplot(tg_env$trtgrps, aes(x = grp, y = est_logodds)) +
    geom_point() +
    scale_y_continuous("Log odds", limits = c(0, max(tg_env$trtgrps$est_logodds)))+
    ggtitle("Log odds")
  
  ggplot(tg_env$trtgrps, aes(x = grp, y = prob_best)) +
    geom_point() +
    scale_y_continuous("Probability", limits = c(0, max(tg_env$trtgrps$prob_best)))+
    ggtitle("Prob best")
  
  ggplot(tg_env$trtgrps, aes(x = grp, y = rand_prob)) +
    geom_point() +
    scale_y_continuous("Probability", limits = c(0, max(tg_env$trtgrps$rand_prob)))+
    ggtitle("Random alloc")
  
  
  
  # only if we have the mcmc samples
  # dfig <- as_tibble(model_mu) %>%
  #   tidyr::gather("trt", "logodds")
  # ggplot(dfig, aes(x = trt, y = logodds))+
  #   geom_jitter(height = 0, width = 0.05)
  #
  # dfig <- as_tibble(model_probs) %>%
  #   tidyr::gather("trt", "prop")
  # ggplot(dfig, aes(x = trt, y = prop))+
  #   geom_jitter(height = 0, width = 0.05)
  
  # hist(model_prop_diffs[, 4])
}




# main loop

interim_seed <- 34432

myseed <- 1
cfg <- get_cfg()

starttime <- Sys.time()
set.seed(myseed)

res <- lapply(1:cfg$nsims, simulate_trial)
dres <- do.call(rbind, res)

endtime <- Sys.time()

difftime(endtime, starttime, units = "hours")

beepr::beep()
w <- warnings()

rdsfilename <- file.path("out",
                         paste0("res-",
                                format(Sys.time(), "%Y-%m-%d-%H-%M-%S"), ".RDS"))

saveRDS(list(results=dres,
             warnings = w,
             starttime = starttime,
             endtime = endtime,
             duration = difftime(endtime, starttime, units = "hours")),
        rdsfilename)
assign("last.warning", NULL, envir = baseenv())

  
  
