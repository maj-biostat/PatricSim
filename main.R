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





pkgs <- c("randomizr", "doParallel","foreach", "dplyr",
          "ggplot2", "rstan", "PatricRCT")

n_sims <- 10
myseed <- 1

starttime <- Sys.time()

# initiate clusters
debug = F
cl <- NA
if(!debug){
  cl <- makeCluster(parallel::detectCores() - 3, outfile="")
  # cl <- makeCluster(3, outfile="")
  registerDoParallel(cl)
} else {
  registerDoSEQ()
}


results <- foreach(i = 1:n_sims,
                   .errorhandling = 'pass',
                   .packages=pkgs
) %dopar%{
  # i = 1
  set.seed(myseed + i)
  res <- simulate_trial(i)
  return(res)
}


dfres <- do.call(rbind, results)
 
endtime <- Sys.time()
difftime(endtime, starttime, units = "hours")

beepr::beep()
w <- warnings()

rdsfilename <- file.path("out", 
                         paste0("res-",
                                format(Sys.time(), "%Y-%m-%d-%H-%M-%S"), ".RDS"))

saveRDS(list(results=dfres,
             warnings = w,
             starttime = starttime,
             endtime = endtime,
             duration = difftime(endtime, starttime, units = "hours")),
        rdsfilename)
assign("last.warning", NULL, envir = baseenv())

if(!debug){
  stopCluster(cl)
}
  
  
