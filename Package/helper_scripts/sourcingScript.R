### Sourcing Files ###

# Loading Function to calculate missing edges, false positive edges, and SHD
source(paste0(path_start,"Desktop/Research/package/test_scripts/local_pc_dist.R"))
# Loading Data Generation Function
source(paste0(path_start,"Desktop/Research/projects/code/data_gen_R.R"))
# Loading function to test local pc for different target sets
source(paste0(path_start,"Desktop/Research/package/helper_scripts/lpc_builds.R"))

# Obtain the names of all the networks
file_vec <- list.files(paste0(ps2,"Desktop/Research/projects/bn_data_generation/networks/rds/"))
nets <- sapply(file_vec, function(x) sub(".rds","",x))

# Function to visualize true DAG in neighborhood against local PC DAG
source(paste0(path_start,"Desktop/Research/package/test_scripts/test_bn.R"))

# Set up result file
if (!dir.exists(paste0(path_start,"Dropbox/Academics/Research/Results/Population")))
  dir.create( paste0(path_start,"Dropbox/Academics/Research/Results/Population") )
