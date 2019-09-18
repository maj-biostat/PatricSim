

##
# In this design we have 4 groups representing:
# 3 days placebo
# 3 days active
# 5 days placebo
# 5 days active

# This is a single bottle design.

# Randomisation is fixed (and balanced) for each arm throughout the trial up to 
# a sample size of 2000.
# No adaptations are in place and a single analysis occurs at n = 2000.
# The response is binary, representing whether a participant recovered
# by day 7 or not. 

# A simple logistic regression model is used for the analysis, thus we 
# do not accommodate censoring. 
# Terms are included for the intercept, the status of the intervention 
# (active vs placebo) and the duration of the intervention (3 and 5 days).

# The final analyses centres on non-inferiority and equivalence.

# library(brms)
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


tg_env <- new.env()

# min samp size for rar by subgroup
rar_min_cell <- 50
non_inf_delta <- 0.1
# want specific ordering.
grp <- factor(c(
  "A5", "A3", "P5", "P3"
), levels = c(
  "A5", "A3", "P5", "P3"
))

A5 <- 1
A3 <- 2
P5 <- 3
P3 <- 4

ACTIVE <- 1
PLACEBO <- 2

# tg_env$trtgrps
# tg_env$rar_active
# tmp2 <- tg_env$trtgrps

interim_seed <- 34432

tg_env$model_code <- rstan::stan_model(file = "logistic_02.stan",
                                       auto_write = TRUE)






scenario <- function(idx = 1){
  # idx = 4
  tg_env$trtgrps <- tibble(
    prob_best = rep(1/length(grp), length(grp)),
    true_mean = rep(0.8, length(grp)),
    prop_rescue = c(0, 0, 0, 0),
    n = 0,
    trt = as.numeric(substr(grp, 1, 1) == "A"),
    grp = grp,
    grplab = as.character(grp),
    # logodds and var used for pbest and rar
    est_logodds = rep(0, length(grp)),
    est_var = rep(0, length(grp))
  )

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
  } else if(idx == 4){
    tg_env$trtgrps$true_mean = rep(0.5, length(grp))
  }else if(idx == 5){
    
    # assume antibiotics are no better than placebo
    # assume 50% have day 7 recovery regardless of the treatment regime
    
    # participants that jump to rescue are considered a failure

    # there will be lower proportion of people recovered at day 3 
    # than day 5 so conclude that there would be more people 
    # jumping to rescue in the 3 day regimes compared to the day 5 
    # regimes
    
    # set 30% of 5 day regime to failure (150)
    # set 40% of 3 day regime to failure (200)
    
    tg_env$trtgrps$true_mean = rep(0.5, length(grp))
    tg_env$trtgrps$prop_rescue = c(0.3, 0.4, 0.3, 0.4)
    
    
  } else if(idx == 6) {
    
    # assume antibiotics do work better than placebo
    # but there is no duration effect 
    # of those on antibiotics, 70% recover by day 7
    # of those on placebo, 50% recover by daty 7
    
    # anyone that jumps to rescue is considered a failure

    # there will be lower proportion of people recovered at day 3 
    # than day 5 so conclude that there would be more people 
    # jumping to rescue in the 3 day regimes compared to the day 5 
    # regimes
    
    tg_env$trtgrps$true_mean = c(0.7, 0.7, 0.5, 0.5)
    tg_env$trtgrps$prop_rescue = c(0, 0, 0.3, 0.4)
    
    
  }
  else if(idx == 7) {

    # assume antibiotics do work better than placebo
    # and there is a duration effect of those on antibiotics
    # Let 70% recover by day 7 if they are on the 5 day regime and
    # 60% recover by day 7 if they are on the 3 day regime.
    
    # For those on placebo assume a 50% recovery by day 7
    
    # there will be lower proportion of people recovered at day 3 
    # than day 5 so conclude that there would be more people 
    # jumping to rescue in the 3 day regimes compared to the day 5 
    # regimes
    
    tg_env$trtgrps$true_mean = c(0.7, 0.6, 0.5, 0.5)
    tg_env$trtgrps$prop_rescue = c(0, 0, 0.3, 0.4)
  } 
  
  # else if(idx == 7) {
  #   # 5 days active better than 3. 
  #   # placebo 3 days reaches for rescue earlier than p5 and therefore
  #   # gets higher proportion recovered.
  #   tg_env$trtgrps$true_mean = c(0.7, 0.6, 0.5, 0.6)
  # } 
}



generate_trial_data <- function() {
  
  
  n <- 2000
  dat <- tibble(
    id = 1:n
  )
  
  # split this up but do not use at the moment.
  # idea is to assess both first 200 as pilot and then all data at 2000.
  grp1 <- randomizr::complete_ra(N = 200, 
                                 num_arms = length(grp),
                                 conditions = grp)
  grp2 <- randomizr::complete_ra(N = n-200, 
                                 num_arms = length(grp),
                                 conditions = grp)
  
  dat$grp <- unlist(list(grp1, grp2))
  
  dat$grp_num <- as.numeric(dat$grp)
  
  message(paste0("Groups produced by randomizr ", length(unique(dat$grp))))

  dat$trt <- as.numeric(substr(dat$grp, 1, 1) == "A")
  dat$dur <- as.numeric(substr(dat$grp, 2, 2) == 5)
  
  # prevents the case where all pts in a given group are all successes
  rngdat <- 1
  grpidx <- as.numeric(dat$grp)
  while(rngdat > 0 & rngdat < 10){
    message(paste0("Data generation attempt ", rngdat))
    dat$y = rbinom(n = n, size = 1, prob = tg_env$trtgrps$true_mean[grpidx])

    cell_n <- dat %>%
      dplyr::group_by(y, grp) %>%
      dplyr::summarise(n = n()) %>%
      dplyr::ungroup() %>%
      dplyr::pull(n)
    
    message(paste0(" ", cell_n))
    
    if(length(cell_n) == length(grp)*2 & all(cell_n > 0)){
      rngdat <- 0
    } else {
      rngdat <- rngdat + 1
    }

    if(rngdat > 10){
      saveRDS(list(dat, cell_n, idx_start, idx_end, grpidx,
                   tg_env$trtgrps$true_mean),
              file = "panic.RDS")
      stop(paste0("Problem generating data."))
    }
  }

  # gmodels::CrossTable(dat$y, dat$grp)
  
  dat$prop_rescue <- tg_env$trtgrps$prop_rescue[dat$grp_num]
  
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
  
  
  # tg_env$df <- generate_trial_data()
  
  #
  # gmodels::CrossTable(tg_env$df$y, tg_env$df$grp)
  
  # participant data
  tmp <- tg_env$df %>%
    dplyr::group_by(trt, dur) %>%
    dplyr::summarise(y = sum(y),
                     n = n()) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(`trt:dur` = trt * dur) %>%
    dplyr::arrange(-trt, -dur)

  Xmat <- as.matrix(tmp %>% dplyr::select(trt, dur, `trt:dur`))
  
  # y = tmp %>% dplyr::pull(y)
  # n = tmp %>% dplyr::pull(n)
  # lm1 <- glm(cbind(y, n-y) ~ trt*dur, data = tmp, family = "binomial"); summary(lm1)
  
  # ## Zmat for estimating mean in each group
  Zmat <- cbind(intercept = 1, Xmat)
  
  model_data <- list(y = tmp %>% dplyr::pull(y),
                     n = tmp %>% dplyr::pull(n),
                     X = Xmat,
                     N = nrow(Xmat), 
                     K = 3,  # excl intercept
                     prior_only = 0)
  
  model_fit <- rstan::sampling(tg_env$model_code, 
                               data = model_data,
                               chains = 1, 
                               iter = 4000,
                               refresh = 4000,
                               seed = interim_seed,
                               control = list(adapt_delta = 0.999),
                               verbose = F)
  
  # model_fit <- rstan::sampling(tg_env$model_code, data = model_data,chains = 1, iter = 1,refresh = 1,seed = interim_seed,control = list(adapt_delta = 0.999),verbose = F)
  
  # print(model_fit)
  
  # log odds of being better by day 7
  model_draws <- as.matrix(model_fit, pars = c("b0", "b"))
  # log odds
  # order is: A5 A3 P5 P3
  model_mu <- model_draws %*% t(Zmat)
  # colMeans(model_mu)

  # because we are looking at non-inferiority and require some form of
  # non-inferiority threshold it makes some sense to convert from log odds
  # to proportions.
  # order is: A5 A3 P5 P3
  model_prop <- apply(model_mu, 2, function(x) exp(x)/(1+exp(x)))
  # sanity check - should be close to the true means used to generate the data
  colMeans(model_prop)

  tg_env$trtgrps$est_logodds <- apply(model_mu, 2, mean)
  tg_env$trtgrps$est_var <- diag(var(model_mu))
  tg_env$trtgrps$prob_best <- p_best(model_mu)
  tg_env$trtgrps$est_prop <- apply(model_prop, 2, mean)
  
  # How about stratifying the analysis by those that jumped to rescue?
  # Or conditioning on duration to jump to rescue (with zero being didn't).
  # Urm....
  # Put this on hold for now.
  
  # Is the short duration active trt non-inferior to long duration trt
  # in terms of the proportion that recover by day 7?
  # There is a worry here that pts on a3 will jump to rescue at day 3 
  # earlier and then continue treatment which would mean that they have had
  # a greater trt duration than those on a5 - unless, of course, those on a5
  # also jumpt to rescue as soon as they finish day 5 :/
  tg_env$mean_a_diff <- mean(model_prop[, A3] - model_prop[, A5])
  tg_env$prob_a3_ni_a5 <-  mean(model_prop[, A3] - model_prop[, A5] > cfg$noninf_threshold)
  tg_env$a3_ni_a5 <- as.numeric(tg_env$prob_a3_ni_a5 > cfg$decision_ni_prob)
  
  # Evaluate the placebos following the same process.
  # Note, to some extent, this should characterise any differential behaviour 
  # with respect to jumping to the rescue medication.
  # Note that p3 could have 4 days of rescue med and p5 could have 2 days of
  # rescue med prior to the day 7 assessment of recovery.
  tg_env$mean_p_diff <- mean(model_prop[, P3] - model_prop[, P5])
  tg_env$prob_p3_ni_p5 <- mean(model_prop[, P3] - model_prop[, P5] > cfg$noninf_threshold)
  tg_env$p3_ni_p5 <- as.numeric(tg_env$prob_p3_ni_p5 > cfg$decision_ni_prob)
  
  # comparison of treatment vs placebo head to head using data pooled across durations
  model_pooled <- matrix(c(model_prop[, A3], model_prop[, A5], 
    model_prop[, P3], model_prop[, P5]), ncol = 2, byrow = F)
  diffs <- model_pooled[, ACTIVE] - model_pooled[, PLACEBO]
  
  tg_env$mean_diff <- mean(model_pooled[, ACTIVE] - model_pooled[, PLACEBO])
  tg_env$prob_a_eq_p <- mean(diffs > -cfg$equiv_threshold & 
                               diffs < cfg$equiv_threshold)
  tg_env$a_eq_p <- as.numeric(tg_env$prob_a_eq_p > cfg$decision_eq_prob)
  
  # Assess equivalence of the difference in differences.
  # In other words is the difference in active (A3 - A5) equivalent to the
  # difference in placebo (P3 - P5)?
  # If there really is no difference between these then we conclude that 
  # providing antibiotic offers no real benefit.
  # To be explicit, in the case where there is no difference between durations 
  # of active treatment we expect that a3 - a5 = 0.
  # Similarly, we anticipate that p3 - p5 = 0 BUT there is also the possibility
  # that p3 - p5 > 0 because, on average, the p3 arm participants would reach 
  # for the rescue medication earlier than the p5 arm. If the rescue medication 
  # is actually helpful then we could presumably see a situation whereby the 
  # proportion that recover by day 7 on p3 is greater than the proportion 
  # that recover by day 7 that were on p5.
  model_prop_diffs <- model_prop[, c(A3, P3)] - model_prop[, c(A5, P5)]
  
  # d1 <- density(model_prop_diffs[,1])
  # d0 <- density(model_prop_diffs[,2])
  # plot(d1, ylim = c(0, max(c(d1$y, d0$y))), xlim = c(-0.2, .1))
  # lines(d0, col = "red")

  # We expect that the difference between the two active arms (A5-A3) will be
  # equivalent to the  difference between the two placebo arms (P5-P3) and so we
  # compute the probability that these two differences are within a threshold of
  # each other. If this probability is high we might be justified in 
  # concluding equivalence.
  diff_in_diffs <- model_prop_diffs[, ACTIVE] - model_prop_diffs[, PLACEBO]
  tg_env$mean_diff_in_diff <- mean(diff_in_diffs)
  tg_env$prob_diff_a_eq_diff_p <- mean(diff_in_diffs > -cfg$equiv_threshold & 
                               diff_in_diffs < cfg$equiv_threshold)
  tg_env$diff_a_eq_diff_p <- as.numeric(tg_env$prob_diff_a_eq_diff_p > cfg$decision_eq_prob)
  
  
  # Alternatively, if this probability was very low we could indicate futile.
  # We could compute the predictive probability of concluding equivalence at the
  # end of the trial.
  

}



simulate_trial <- function(id_trial = 1){
  
  # id_trial = 1
  
  message("")
  message("")
  
  message(paste0("###########################################"))
  message(paste0("  TRIAL ", id_trial, " of ", cfg$nsims, " SCENARIO ", cfg$scenarioid))
  message(paste0("###########################################"))
  
  # scenario(cfg$scenarioid)
  scenario(cfg$scenarioid)

  # reset data
  tg_env$df <- generate_trial_data()
  
  # gmodels::CrossTable(tg_env$df$rescue_sought, tg_env$df$grp)
  # gmodels::CrossTable(tg_env$df$y, tg_env$df$grp)
  
  tg_env$trtgrps$n <- tg_env$df %>%
    dplyr::group_by(grp) %>%
    dplyr::summarise(n = n()) %>%
    dplyr::ungroup() %>%
    # if one of the groups has not been randomised any patients
    # then explicitly set it to zero rather than having it missing
    dplyr::right_join(tibble(grp = tg_env$trtgrps$grp), by = "grp") %>%
    dplyr::mutate(n = ifelse(is.na(n), 0, n)) %>%
    dplyr::pull(n)
  
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

  
  
