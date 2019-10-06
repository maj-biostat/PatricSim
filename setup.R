

get_cmdline_opts <- function(){
  
  option_list <- list(
    optparse::make_option(c("-s", "--scenario"), 
                          type = "integer", default = NULL,
                          help = "scenario", metavar = "integer"),
    optparse::make_option(c("-n", "--nsims"), 
                          type = "integer", default = NULL,
                          help = "number of simulations", 
                          metavar = "integer")
  )
  
  opt_parser <- OptionParser(option_list = option_list);
  opt <- parse_args(opt_parser);
  opt
  
}


# Strictly speaking I do not require the config file but will leave in 
# for the time being
get_cfg <- function(cfgfile = "cfg.yaml"){
  
  opt <- get_cmdline_opts()
  print(opt)

  tt <- tryCatch(configtmp <- read.config(file = file.path(getwd(), cfgfile)),
                 error=function(e) e,
                 warning=function(w) w)
  ifelse(is(tt,"warning") | is(tt,"error"),"Configuration Warning/Error.
         Please ensure configuration file has terminating empty line.",
         "Configuration File Loaded OK")
  
  l <- list()
  
  if("scenario" %in% names(opt)){
    l$scenarioid <- opt$scenario
  } else {
    l$scenarioid <- tt$scenarioid
  }
  
  if("nsims" %in% names(opt)){
    l$nsims <- opt$nsims
  } else {
    l$nsims <- tt$nsims
  }

  
  if("dose" %in% names(opt)){
    l$dose <- 0
  } else {
    
    l$dose <- tt$dose
    
  }
  
  
  l$n_per_trt <- tt$n_per_trt
  l$p_ni_thresh <- tt$p_ni_thresh
  l$p_ni_deicsion_thresh <- tt$p_ni_deicsion_thresh
  l$p_sup_deicsion_thresh <- tt$p_sup_deicsion_thresh
  
  l
}

  