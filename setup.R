

get_cmdline_opts <- function(){
  
  option_list <- list(
    optparse::make_option(c("-s", "--scenario"), 
                          type = "integer", default = NULL,
                          help = "scenario", metavar = "integer")
  )
  
  opt_parser <- OptionParser(option_list = option_list);
  opt <- parse_args(opt_parser);
  opt
  
}


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


  l
}

  