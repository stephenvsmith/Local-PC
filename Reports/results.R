##################################################################

# Author: Stephen Smith
# Date: 1/13/20
# Description: This file is for compiling results from 
# the simulations

##################################################################

library(ggplot2)

##### First, we grab results with tolerance = 0.005

setwd("~/Desktop/Research/Results/saved/sample/tol=0.005")
files <- list.files()

res <- list()
for (f in files){
  res[[f]] <- read.table(f)
}

all_results <- do.call(rbind,res)

p <- ggplot(data = all_results, aes(x = Time_to_Build))+geom_histogram()
p + facet_wrap(~Network,scales = "free")

p <- ggplot(data = all_results, aes(x = Undirected.Missing))+geom_histogram()
p + facet_wrap(~Network)

p <- ggplot(data = all_results, aes(x = Directed.Missing))+geom_histogram()
p + facet_wrap(~Network)

p <- ggplot(data = all_results, aes(x = Total.Missing))+geom_histogram()
p + facet_wrap(~Network)

p <- ggplot(data = all_results, aes(x = Undirected.Added))+geom_histogram()
p + facet_wrap(~Network)

p <- ggplot(data = all_results, aes(x = Directed.Added))+geom_histogram()
p + facet_wrap(~Network)

p <- ggplot(data = all_results, aes(x = Total.Added))+geom_histogram()
p + facet_wrap(~Network)

p <- ggplot(data = all_results, aes(x = Wrong_Direction))+geom_histogram()
p + facet_wrap(~Network)

p <- ggplot(data = all_results, aes(x = Directed_Undirected))+geom_histogram()
p + facet_wrap(~Network)

p <- ggplot(data = all_results, aes(x = Undirected_Directed))+geom_histogram(bins=10)
p + facet_wrap(~Network)


##############################################################################################33

# All Results with tol = 0.01

setwd("~/Desktop/Research/Results/sample")

files <- list.files()

res2 <- list()
for (net in files){
  if (dir.exists(paste0("./",net))){
    if (file.exists(paste0(net,"/",net," results.txt"))) {
      res2[[net]] <- read.table(paste0(net,"/",net," results.txt"))
    }
  }
}

all_results2 <- do.call(rbind,res2)

  
