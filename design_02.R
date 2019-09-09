
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

# define treatment groups, initial allocation probability and true means

tg_env <- new.env()
tg_env$trtgrps <- tibble(
  prob_best = rep(1/2, 2),
  rand_prob = rep(1/2, 2),
  true_mean = rep(0.8, 2),
  n = 0,
  grp = as.factor(paste0("T", 1:2)),
  grplab = as.character(paste0(c(7, 0), "-day")),
  # logodds and var used for pbest and rar
  est_logodds = rep(0, 2),
  est_var = rep(0, 2),
  est_prop = rep(0, 2),
  prob_ni = rep(0, 2),
  trt_ni = rep(0, 2)
)
tg_env$rar_active <- 0

rownames(tg_env$trtgrps) <- NULL







generate_trial_data <- function(idx_start = 1, idx_end = 50) {
  
  
  n <- idx_end - idx_start + 1
  dat <- tibble(
    id = idx_start:idx_end
  )
  
  n_grps <- max(as.numeric(tg_env$trtgrps$grp))
  
  if(!tg_env$rar_active){
    dat$grp <- randomizr::complete_ra(N = n, num_arms = 2)
  } else {
    dat$grp <- randomizr::complete_ra(N = n,
                                      num_arms = n_grps,
                                      prob_each = tg_env$trtgrps$rand_prob)
  }
  
  # prevents the case where all pts in a given group are all successes
  rngdat <- T
  its <- 1
  while(rngdat){
    dat$y = rbinom(n = n, size = 1, prob = tg_env$trtgrps$true_mean[dat$grp])
    
    cell_n <- dat %>%
      dplyr::group_by(y, grp) %>%
      dplyr::summarise(n = n()) %>%
      dplyr::ungroup() %>%
      dplyr::pull(n)
    
    if(length(cell_n) == n_grps*2 & all(cell_n > 0)){
      rngdat <- F
    }
    
    its <- its + 1
    
    if(its > 10){
      dat
      table(dat$y, dat$grp)
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



fit_stan <- function(df){
  
  # participant data
  Xmat <- model.matrix(y ~ 0 + grp, data = df)
  
  n_grps <-  max(as.numeric(tg_env$trtgrps$grp))
  
  model_code <- rstan::stan_model(file = "logistic.stan", auto_write = TRUE)
  model_data <- list(y = df$y,
                     X = Xmat,
                     N = length(df$y),
                     K = n_grps,
                     prior_only = 0)
  
  model_fit <- rstan::sampling(model_code, data = model_data,
                               chains = 1, iter = 2000,
                               seed = 652732,
                               control = list(adapt_delta = 0.999))
  
  # log odds of being better by day 7
  model_mu <- as.matrix(model_fit, pars = c("b"))
  model_props <- apply(model_mu, 2, function(x) exp(x)/(1+exp(x)))
  colnames(model_props) <- paste0("b", 1:n_grps)
  tg_env$trtgrps$est_prop <- colMeans(model_props)
  
  # comparisons
  model_prop_diffs <- do.call(cbind, lapply(1:n_grps, function(x)
    model_props[,1] - model_props[,x]))
  
  tg_env$trtgrps$est_logodds <- apply(model_mu, 2, mean)
  tg_env$trtgrps$est_var <- diag(var(model_mu))
  tg_env$trtgrps$prob_best <- p_best(model_mu)
  tg_env$trtgrps$prob_ni <- colMeans(model_prop_diffs < non_inf_delta)
  
  tg_env$trtgrps$rand_prob <- sqrt((tg_env$trtgrps$prob_best *
                                      tg_env$trtgrps$est_var)/(tg_env$trtgrps$n + 1))
  
  tg_env$trtgrps$rand_prob <- tg_env$trtgrps$rand_prob / sum(tg_env$trtgrps$rand_prob)
  
}





simulate_trial <- function(id_trial){
  
  scenario(1)
  
  looks <- seq(100,500,by=50)
  
  df <- tibble()
  
  for(idx_intrm in 1:length(looks)){
    
    # idx_intrm=1
    # idx_intrm=idx_intrm+1
    message(paste0("Interim ", idx_intrm))
    df <- rbind(df,
                generate_trial_data(idx_start = nrow(df)+1,
                                    idx_end = looks[idx_intrm]))
    
    # rar enabled from the second interim
    tg_env$rar_active <- 1
    
    # cumulative total of participants assigned to each group
    tg_env$trtgrps$n <- df %>%
      dplyr::group_by(grp) %>%
      dplyr::summarise(n = n()) %>%
      dplyr::ungroup() %>%
      # if one of the groups has not been randomised any patients
      # then explicitly set it to zero rather than having it missing
      dplyr::right_join(tibble(grp = tg_env$trtgrps$grp), by = "grp") %>%
      dplyr::mutate(n = ifelse(is.na(n), 0, n)) %>%
      dplyr::pull(n)
    
    # table(df$subX)
    
    fit_stan(df)
    
    trts_ni <- c(tg_env$trtgrps$prob_ni > tg_env$ni_thresh)
    if(all(trts_ni)){
      tg_env$trtgrps$trt_ni <- trts_ni
      tg_env$trtgrps$trt_ni[1] <- NA
      message("All treatments NI", paste0(trts_ni, sep= " "))
      break
    }
    
  }
  
  
  cbind(id = id_trial, tg_env$trtgrps)
  
  
}




















