source("fun_graph.r")
source("fun_network.r")

##----- input data:  a network matrix and a list of p values from association tests; p values should be in the same order as genes in the network.csv file
net <- read.table("network.csv",sep=",",header=T, row.names=1)
pval <- read.table("pval.txt", header=F)

##----- set priors and MCMC parameters
priorPar1 <- c(-2, 0.20, 0.05)
priorPar2 <- c(3,1,10,1)
mcmcNum <- 100
mcmcBurnin <- 20
mstart <- 30

##--------- start preprocessing data --------
ngene <- ncol(net)
genesymbol <- colnames(net)
edges <- expand.grid(genesymbol,genesymbol,stringsAsFactors =F)[as.vector(net)==1,]
colnames(edges) <- c("genesymbol1","genesymbol2")
netwk <- procGraph(data.frame(genesymbol=genesymbol),edges)
trueEdges <- edgeClass(netwk$trueEdges,0.5)
nNodes <- attr(trueEdges,"nNodes")
nEdges <- nrow(trueEdges)
nodedata <- netwk$nodes
hist(nodedata$nNei,breaks=10)
obsAsso <- pval.class(pval[,1])

##----- run Gibbs sampler ----
gsrst <- Gibbs.mstart( mstart, mcmcNum,  trueEdges, priorPara=priorPar1, priorPara2=priorPar2, nodesData=obsAsso,
                      initLabel=NULL,sampleFrmPrior=F)
##---- analyse MCMC data ----
ana <- gmstart.rst(trueEdges,gsrst,mcmcBurnin)

#---- results from the MCMC ----
## MAP
map <- ana$mostFreq
##Posterior mean
postMean <- ana$postMean
##MMP
mmp <- (postMean>=0.5)*2-1
##MCP
mcp <- icm(nodesData=obsAsso, trueEdges, priorPara=priorPar1, priorPara2=priorPar2,initLabel=NULL,
           searchMethod=c("DFS","BFS","random"), printI=T) 

