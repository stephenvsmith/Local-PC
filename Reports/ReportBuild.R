res <- read.table("~/Desktop/Research/Results/sample/alarm/alarm results.txt")
res

summary(res$Same_V_Structures)

setwd("~/Desktop/Research/Results/sample")
dirs <- list.files()
for (d in dirs){
  if (file.exists(paste0("./",d,"/",d," results.txt"))){
    #dir.create(paste0("~/Desktop/Research/Results/saved/sample/",d))
    file.copy(paste0("./",d,"/",d," results.txt"),paste0("~/Desktop/Research/Results/saved/sample/",d," results.txt"))
  }
}
