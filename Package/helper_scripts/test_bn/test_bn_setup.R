#test_bn_setup

file <- paste0(path_start,"Desktop/Research/projects/bn_data_generation/networks/rds/",net,".rds")
names <- names(readRDS(file))
t_dag_original <- empty.graph(names)
amat(t_dag_original) <- true_dag

cat("Network: ",net,"\n",
    file = paste0(path_start,"Dropbox/Academics/Research/Build Notes/test_bn_notes.txt"),append = FALSE)
result <- local_pc2(true_dag = true_dag,target = target,verbose = FALSE,verbose_small = FALSE,pop = TRUE)

cat("Local PC Completed Running\n",
    file = paste0(path_start,"Dropbox/Academics/Research/Build Notes/test_bn_notes.txt"),append = TRUE)
