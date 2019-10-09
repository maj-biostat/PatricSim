






load_results_data <- function(dirname, p_sup = 0.9, p_ni = 0.6){
  
  # dirname = "out05c"
  # dirname = "out05e"
  
  dsup <- data.frame()
  dni <- data.frame()
  
  drecov <- data.frame()
  drecov_lwr <- data.frame()
  drecov_upr <- data.frame()
  dbias <- data.frame()
  
  dtrue_mu <- data.frame()
  dcompliance <- data.frame()
  
  fout <- list.files(dirname)
  
  for(i in 1:length(fout)){
    l <- readRDS(paste0(dirname, "/", fout[i]))
    
    # superiority data
    
    # test whether max is superior to all other doses
    dtmp <- cbind(l$dsup, apply(l$dsup[, 1:ncol(l$dsup)], 2, function(x){ x > p_sup}))
    
    n_doses <- length(l$trtgrps$trtgrps$dose)
    doses <- l$trtgrps$trtgrps$dose
    min_dose <- min(l$trtgrps$trtgrps$dose)
    max_dose <- max(l$trtgrps$trtgrps$dose)
    mu_lwr <- min(l$trtgrps$trtgrps$true_mu)
    mu_upr <- max(l$trtgrps$trtgrps$true_mu)
    compliance <- l$trtgrps$trtgrps$compliance
    
    
    dcompliance <- rbind(dcompliance, compliance)
    
    
    
    tmp <- c(scenario = l$scenario, 
             range = l$trtgrps$p_range,
             mid = l$trtgrps$location,
             mu_lwr = mu_lwr,
             mu_upr = mu_upr,
             n_doses = n_doses,
             min_dose = min_dose,
             max_dose = max_dose,
             
             n_per_trt = l$n_per_trt,
             colMeans(dtmp[, n_doses:ncol(dtmp)]))
    
    # eta columns are the proportion of times that the proportion that recovered on max
    # duration was greater than the proportion recovered on the 1 day, 2 days, 3 days etc
    dsup <- rbind(dsup, tmp)
    
    names(dsup) <- names(tmp)
    names(dsup)[(ncol(dsup)-n_doses+2):ncol(dsup)] <- paste0("D", doses[-length(doses)])
    
    # non-inferiority
    # dtmp <- cbind(l$dni, apply(l$dni[, 1:ncol(l$dni)], 2, mean))
    
    dtmp <- cbind(l$dni, apply(l$dni[, 1:ncol(l$dni)], 2, function(x){ x > p_ni}))
    
    tmp <- c(scenario = l$scenario, 
             range = l$trtgrps$p_range,
             mid = l$trtgrps$location,
             mu_lwr = mu_lwr,
             mu_upr = mu_upr,
             n_doses = n_doses,
             min_dose = min_dose,
             max_dose = max_dose,
             n_per_trt = l$n_per_trt,
             colMeans(dtmp[, n_doses:ncol(dtmp)]))
    
    # eta columns are the proportion of times that the proportion that recovered on max
    # duration was non-inferior to the proportion recovered on the 1 day, 2 days, 3 days etc
    dni <- rbind(dni, tmp)
    
    names(dni) <- names(tmp)
    names(dni)[(ncol(dni)-n_doses+2):ncol(dni)] <- paste0("D", doses[-length(doses)])
    
    # prob recovery
    tmp <- c(scenario = l$scenario, 
             range = l$trtgrps$p_range,
             mid = l$trtgrps$location,
             mu_lwr = mu_lwr,
             mu_upr = mu_upr,
             n_doses = n_doses,
             min_dose = min_dose,
             max_dose = max_dose,
             n_per_trt = l$n_per_trt,
             colMeans(l$dprob))
    
    drecov <- rbind(drecov, tmp)
    names(drecov) <- names(tmp)
    names(drecov)[(ncol(drecov)-n_doses+1):ncol(drecov)] <- paste0("D", doses)
    
    tmp <- c(scenario = l$scenario, 
             range = l$trtgrps$p_range,
             mid = l$trtgrps$location,
             mu_lwr = mu_lwr,
             mu_upr = mu_upr,
             n_doses = n_doses,
             min_dose = min_dose,
             max_dose = max_dose,
             n_per_trt = l$n_per_trt,
             colMeans(l$dprob_lwr))
    
    drecov_lwr <- rbind(drecov_lwr, tmp)
    names(drecov_lwr) <- names(tmp)
    names(drecov_lwr)[(ncol(drecov_lwr)-n_doses+1):ncol(drecov_lwr)] <- paste0("D", doses)
    
    tmp <- c(scenario = l$scenario, 
             range = l$trtgrps$p_range,
             mid = l$trtgrps$location,
             mu_lwr = mu_lwr,
             mu_upr = mu_upr,
             n_doses = n_doses,
             min_dose = min_dose,
             max_dose = max_dose,
             n_per_trt = l$n_per_trt,
             colMeans(l$dprob_upr))
    
    drecov_upr <- rbind(drecov_upr, tmp)
    names(drecov_upr) <- names(tmp)
    names(drecov_upr)[(ncol(drecov_upr)-n_doses+1):ncol(drecov_upr)] <- paste0("D", doses)
    # prob recovery end
    
    tmp <- cbind(scenario = l$scenario, 
                 range = l$trtgrps$p_range,
                 mid = l$trtgrps$location,
                 mu_lwr = mu_lwr,
                 mu_upr = mu_upr,
                 n_doses = n_doses,
                 min_dose = min_dose,
                 max_dose = max_dose,
                 dose = l$trtgrps$trtgrps$dose,
                 true_mu = l$trtgrps$trtgrps$true_mu
    )
    
    dtrue_mu <- rbind(dtrue_mu, tmp)
    
    
    # bias
    tmp <- c(scenario = l$scenario, 
             range = l$trtgrps$p_range,
             mid = l$trtgrps$location,
             mu_lwr = mu_lwr,
             mu_upr = mu_upr,
             n_doses = n_doses,
             min_dose = min_dose,
             max_dose = max_dose,
             n_per_trt = l$n_per_trt,
             colMeans(l$dbias))
    
    dbias <- rbind(dbias, tmp)
    names(dbias) <- names(tmp)
    names(dbias)[(ncol(dbias)-n_doses+1):ncol(dbias)] <- paste0("D", doses)
    
  }
  
  rownames(dsup) <- NULL
  rownames(dni) <- NULL
  
  names(dcompliance) <- paste0("Dose ", doses)
  
  dcompliance <- cbind(scenario = 1:nrow(dcompliance), dcompliance)
  
  dsup <- dsup %>% dplyr::arrange(scenario)
  dni <- dni %>% dplyr::arrange(scenario)
  drecov <- drecov %>% dplyr::arrange(scenario)
  drecov_lwr <- drecov_lwr %>% dplyr::arrange(scenario)
  drecov_upr <- drecov_upr %>% dplyr::arrange(scenario)
  dbias <- dbias %>% dplyr::arrange(scenario)
  dtrue_mu <- dtrue_mu %>% dplyr::arrange(scenario)
  
  
  return(list(dsup = dsup,
              dni = dni, 
              drecov = drecov,
              drecov_lwr = drecov_lwr,
              drecov_upr = drecov_upr,
              dbias = dbias,
              dtrue_mu = dtrue_mu,
              
              n_doses = n_doses,
              doses = doses,
              min_dose = min_dose,
              max_dose = max_dose,
              dcompliance = dcompliance))
  
  
}


