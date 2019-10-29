### Save DAG photos

# Create a File for the Network
if (!file.exists(paste0("./",net))) dir.create(paste0("./",net))
setwd(paste0("./",net))

# Save File of True DAG
file <- paste0(ps2,"Desktop/Research/projects/bn_data_generation/networks/rds/",net,".rds")
names <- names(readRDS(file))
t_dag <- empty.graph(names)
rownames(true_dag) <- names
colnames(true_dag) <- names
amat(t_dag) <- true_dag
png(file = "true_dag.png")
par(mfrow=c(1,1))
bnlearn::graphviz.plot(t_dag,main = "DAG")
dev.off()

# Save File of True CPDAG
tcp_dag <- cpdag(t_dag)
png(file = "true_cp_dag.png")
par(mfrow=c(1,1))
bnlearn::graphviz.plot(t_dag,main = "CPDAG")
dev.off()

