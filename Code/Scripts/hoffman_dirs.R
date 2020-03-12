### Directories for Hoffman

hoffman_dirs <- function(vars){
  # Define paths for storing
  vars$home_dir <- '/u/scratch/s/stephens'
  vars$package_loc <- paste0(vars$home_dir,'/Hoffman/Packages/LocalPC')
  vars$data_gen_file <- paste0(vars$home_dir,'/Hoffman/data_gen.R')
  vars$result_dir <- paste0(vars$home_dir,'/Hoffman/Results/',pop_string,'/')
  
  if (!vars$pop) {
    vars$result_dirs <- sapply(vars$alpha_vals,function(a){
      result_dir <- paste0(vars$result_dir,'alpha=',a)
      if (!dir.exists(vars$result_dir)){
        dir.create(vars$result_dir)
      }
      return(result_dir)
    })
    names(vars$result_dirs) <- paste0("alpha=",vars$alpha_vals)
  }
  vars$simulation_dir <- paste0(vars$home_dir,'/Hoffman/Simulations/')
  vars$rds_dir <- paste0(vars$home_dir,'/Hoffman/rds')
  
  return(vars)
}