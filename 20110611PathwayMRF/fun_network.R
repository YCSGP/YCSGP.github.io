truncCondProb <- F
truncCondProbRange <- c(0.01,0.99)

pval.class <- function(pv) {
  pv[pv==0] <- .Machine$double.eps
  nodeD <- pv
  class(nodeD) <- "pval"
  attr(nodeD,'p') <- pv
  attr(nodeD,'y') <-  qnorm(1-pv/2)
  nodeD
}

"[.norm" <- function(x,i) {
  xx <- .subset(x,i);
  aa <- attributes(x);
  aa$names <- .subset(aa$names,i);
  aa$p <- .subset(aa$p,i);
  attributes(xx) <- aa;
  xx
}

"[.pval" <- function(x,i) {
  xx <- .subset(x,i);
  aa <- attributes(x);
  aa$names <- .subset(aa$names,i);
  aa$p <- .subset(aa$p,i);
  aa$y <- .subset(aa$y,i);
  attributes(xx) <- aa;
  xx
}

edgeClass <- function(edges,wpower,nNodes=NULL) {
  if(is.null(nNodes)) {
    nNodes <- length(union(edges[,1],edges[,2]))
  }
 
  neighbors <- allNeighbour(nNodes,edges)
  w <- (sapply(neighbors,length))^wpower
  
  attr(edges,'wpower') <- wpower
  attr(edges,'neighbors') <- neighbors
  attr(edges,'nNodes') <- nNodes
  attr(edges,'w') <- w
  edges
}



likpw <- function(nodesData,...) {
  UseMethod("likpw", nodesData)
}

likpw.norm <- function(nodesData, nodesLabel) {
  dnorm(as.numeric(nodesData), nodesLabel==1, attr(nodesData,"std"))
}


likpw.pval <- function( nodesData, nodesLabel, priorPara2=NULL) {
  y <- attr(nodesData,'y')
  mubar <- priorPara2[1]
  a <- priorPara2[2]
  nu <- priorPara2[3]
  lambda <- priorPara2[4]
  scale0 <- sqrt(lambda*(a+1)/a)
  ifelse(nodesLabel==1, dt((y-mubar)/scale0,nu)/scale0  , dnorm(y,0,1) )
}


estPriorPara2 <- function(nodesData,...) {
  UseMethod("estPriorPara2", nodesData)
}

estPriorPara2.norm <- function(nodesData) {
  return(0);
}


estPriorPara2.pval <- function( nodesData) {
  y <- attr(nodesData,'y')
  n <- length(y)
  ysort <- sort(y,decreasing=T)
  sampH1 <- ysort[1:max(n*0.2,1)] 
  sampvar <- var(y)
  H1var <- var(sampH1)
  if(is.na(H1var)) H1var <- sampvar
  return(c( mubar=max(3,mean(sampH1)), a= max(1,sampvar/H1var), nu=10, lambda=max(1,H1var)))
}


label2num <- function(label,nNodes=length(label)) {
  bases <- 2^((nNodes-1):0)
  sum((label==1)*bases)
}
num2label <- function(num,nNodes) {
  bases <- 2^((nNodes-1):0)
  label <- rep(0,nNodes)
  for ( kk in 1:nNodes) {
    label[kk] <- num %/% bases[kk]
    num <- num %% bases[kk]
  }
  2*label-1
}






graphPrior <- function(nodesLabel,  edges, h1, tau00, tau11) {
  
  w <- attr(edges,'w')
  edge00 <- nodesLabel[edges[,1]] ==-1 & nodesLabel[edges[,2]]==-1
  edge11 <- nodesLabel[edges[,1]] == 1 & nodesLabel[edges[,2]]== 1
  U <- -tau00*sum( w[ edges[edge00,1] ] + w[ edges[edge00,2] ]) - tau11*sum( w[ edges[edge11,1] ] + w[ edges[edge11,2] ])- h1*sum(nodesLabel==1)
  exp(-U)
}


enumeratePri <- function(trueEdges, h1vec, tau00vec, tau11vec=NULL){
  nodesid <- sort(union(trueEdges[,1],trueEdges[,2]))
  nNodes <- length(nodesid)
  tau.eql <- is.null(tau11vec) 
  if(tau.eql) {
    totalrows <- length(h1vec)*length(tau00vec)
    paramat <- as.matrix(expand.grid(h1vec, tau00vec))
    paramat <- cbind(paramat, paramat[,2])
    colnames(paramat) <- c("h1","tau00","tau11")
  }
  else {
    totalrows <- length(h1vec)*length(tau00vec)*length(tau11vec)
    paramat <- as.matrix(expand.grid(h1vec, tau00vec, tau11vec))
    colnames(paramat) <- c("h1","tau00","tau11")
  }

  networkLab <- matrix(0, nrow=nNodes, ncol=2^nNodes, dimnames=list(paste("node",1:nNodes,sep=""),paste("config",1:2^nNodes,sep=""))) 
  networkPri <- matrix(0, nrow=totalrows, ncol= 2^nNodes, dimnames=list(NULL,paste("config",1:2^nNodes,sep="")) )

  
  for(i in 1:(2^nNodes)) { 
    networkLab[,i] <- num2label(i-1,nNodes)
  }
  
  for (rowi in 1:totalrows) {
    
    h1 <- paramat[rowi,"h1"]
    tau00 <- paramat[rowi, "tau00"]
    tau11 <- paramat[rowi, "tau11"]
    for(i in 1:(2^nNodes)) { 
      networkPri[rowi,i] <- graphPrior(networkLab[,i],trueEdges, h1, tau00, tau11)
      networkPri[rowi,i] <- graphPrior(networkLab[,i],trueEdges, h1, tau00, tau11)
    }
  }
  
  normConst <- apply(networkPri,1,sum)
  networkPriUnnorm <- networkPri
  networkPri <- sweep(networkPri,1,apply(networkPri,1,sum),FUN="/")
  
  nodemarg <- t(apply(networkPri, 1, function(prob) apply(networkLab,1,function(nodei) sum((nodei==1)*prob))))
  colnames(nodemarg) <- paste("node",1:nNodes,sep="")
  list(networkLab=networkLab, paramat=paramat, networkPri=networkPri, networkPriUnnorm=networkPriUnnorm, normConst=normConst, nodemarg=nodemarg)
}


enumeratePost <- function(trueEdges, nodesData, h1vec, tau00vec, tau11vec=NULL, priorPara2){
  nodesid <- sort(union(trueEdges[,1],trueEdges[,2]))
  nNodes <- length(nodesid)
  tau.eql <- is.null(tau11vec) 
  if(tau.eql) {
    totalrows <- length(h1vec)*length(tau00vec)
    paramat <- as.matrix(expand.grid(h1vec, tau00vec))
    paramat <- cbind(paramat, paramat[,2])
    colnames(paramat) <- c("h1","tau00","tau11")
  }
  else {
    totalrows <- length(h1vec)*length(tau00vec)*length(tau11vec)
    paramat <- as.matrix(expand.grid(h1vec, tau00vec, tau11vec))
    colnames(paramat) <- c("h1","tau00","tau11")
  }

  networkLab <- matrix(0, nrow=nNodes, ncol=2^nNodes, dimnames=list(paste("node",1:nNodes,sep=""),paste("config",1:2^nNodes,sep=""))) 
  networkPost <- matrix(0, nrow=totalrows, ncol= 2^nNodes, dimnames=list(NULL,paste("config",1:2^nNodes,sep="")) )

  
  for(i in 1:(2^nNodes)) { 
    networkLab[,i] <- num2label(i-1,nNodes)
  }

  lik.1 <- likpw(nodesData, rep(1,nNodes), priorPara2)
  lik.0 <- likpw(nodesData, rep(-1,nNodes))
  

  for (rowi in 1:totalrows) {
    
    h1 <- paramat[rowi,"h1"]
    tau00 <- paramat[rowi, "tau00"]
    tau11 <- paramat[rowi, "tau11"]
    for(i in 1:(2^nNodes)) { 
      
    networkPost[rowi,i] <- prod(ifelse(networkLab[,i]==1, lik.1, lik.0)) * graphPrior(networkLab[,i],trueEdges, h1, tau00, tau11)
    }
  }
  
  normConst <- apply(networkPost,1,sum)
  networkPostUnnorm <- networkPost
  networkPost <- sweep(networkPost,1,apply(networkPost,1,sum),FUN="/")
  
  nodemarg <- t(apply(networkPost, 1, function(prob) apply(networkLab,1,function(nodei) sum((nodei==1)*prob))))
  colnames(nodemarg) <- paste("node",1:nNodes,sep="")
  list(networkLab=networkLab, paramat=paramat, networkPost=networkPost, networkPostUnnorm=networkPostUnnorm, normConst=normConst, nodemarg=nodemarg)
}



simulateData <-  function(classnm,nodesMean, nodesStd,alpha=0.05,m=NULL,beta=NULL){
  if (classnm=="pval") {
    nodeD <- simulatePval(nodesMean,alpha,m,beta)
  }
  else {
    nodeD <- rnorm(length(nodesMean),nodesMean, nodesStd)
    pvalue <- 2*(1-pnorm(abs(nodeD- (-1) ), 0, nodesStd)) 
    
    class(nodeD) <- classnm
    attr(nodeD,'std') <- nodesStd
    attr(nodeD,'p') <- pvalue
    
  }
  nodeD
}

simulatePval <-  function(nodesLabel, mu,sigma) {
  
  pvalue <- sapply(nodesLabel,function(lbl) if(lbl==1) {rnm <- rnorm(1,mu,sigma); if(rnm>0) 2*(1-pnorm(rnm,0,1)) else 2*pnorm(rnm,0,1) } else runif(1))
  pvalue[pvalue==0] <- .Machine$double.eps
  nodeD <- pvalue
  class(nodeD) <- "pval"
  attr(nodeD,'p') <- pvalue
  attr(nodeD,'y') <-  qnorm(1-pvalue/2) 
  nodeD
}






Gibbs.t <- function( iteNum, edges, priorPara, priorPara2, temperature, initLabel=NULL, nodesData=NULL, sampleFrmPrior=F){ 
  if(!sampleFrmPrior && is.null(nodesData)) stop("Observed data for each node is required!")
  
  
  
  MRF.condProb <- function(nodei,nodesLabel) { 
    i.neighbor <- neighbors[[nodei]] 
    i.neighbor.0 <- i.neighbor[nodesLabel[i.neighbor]==-1]
    i.neighbor.1 <- i.neighbor[nodesLabel[i.neighbor]== 1]    
    logitp <- ( h1 + tau11*( w[nodei]*sum(nodesLabel[i.neighbor]== 1)+sum(w[ i.neighbor.1])) -
                     tau00*( w[nodei]*sum(nodesLabel[i.neighbor]==-1)+sum(w[ i.neighbor.0]))   )/temperature
    1/(1+exp(-logitp))
  }
  
  
  Gibbs.onestep.sampleFrmPrior <- function() {
    currentPos <<- currentPos+1
    nextNodesLabel <- currentNodesLabel    
    for (nodei in sample(nNodes)) { 
      nextNodesLabel[nodei] <- 2*rbinom(1,1, MRF.condProb(nodei, nextNodesLabel))-1
    }
    currentNodesLabel <<- nextNodesLabel
    allPostLik[currentPos] <<- graphPrior(currentNodesLabel,edges,h1, tau00, tau11) 
    allNodesLabel[currentPos,] <<- currentNodesLabel
    allPostLik[currentPos]
  }
  
  Gibbs.onestep <- function() {
    currentPos <<- currentPos+1
    nextNodesLabel <- currentNodesLabel
    for (nodei in sample(nNodes) ) { 
      condProb <- MRF.condProb(nodei, nextNodesLabel)
      if(truncCondProb) condProb <- min(max(condProb,truncCondProbRange[1]),truncCondProbRange[2]) 
      prob <- (likRatio[nodei]*condProb) / (likRatio[nodei]*condProb + (1-condProb))
      nextNodesLabel[nodei] <- 2*rbinom(1,1, prob)-1
    }
    currentNodesLabel <<- nextNodesLabel
    allPostLik[currentPos] <<- graphPrior(currentNodesLabel,edges,h1, tau00, tau11)*prod( ifelse( currentNodesLabel==1,lik.1,lik.0 ))
    allNodesLabel[currentPos,] <<- currentNodesLabel
    allPostLik[currentPos]
  }
  getChainLast <- function() {
    list(nodeLabel=allNodesLabel[currentPos,], postLik=allPostLik[currentPos]);
  }   
    setChainLast <- function(valList) {
      allNodesLabel[currentPos,] <<- valList$nodeLabel
    allPostLik[currentPos] <<- valList$postLik
    return(NULL);
  } 
  
  getTrace <- function(outfile="", header=T,append=F){
    if( !append) cat(file=outfile);
    if (header)
      cat("Temperature","model_idx", dataObj$xNames, "xnum", "freq","lambda","pMarg",
          "post", "\n", file=outfile,append=T);

    for ( i in 1:currentPos) {
      gammaM <- gammaN2M(gammaChain[i]);
      cat(temperature, gammaChain[i], gammaM, sum(gammaM),freqChain[i], lambdaChain[i],datMargPChain[i],
          postChain[i], "\n", file=outfile,append=T)
    }
    
    length(gammaChain) <<- length(lambdaChain) <<- length(phiChain) <<- 
      length(datMargPChain) <<- length(postChain) <<- length(freqChain) <<- currentPos; 
    save( currentPos, temperature, gammaChain, lambdaChain, phiChain,
         datMargPChain, postChain, freqChain,
         file=paste(outfile,"T", temperature,".rdata", sep=""), compress=T);
  }
  
  getMostVisited <- function(outfile="", keepN=currentPos, header=T,append=F){
    if( !append) cat(file=outfile);
    cat("Accept rate:", acceptRate/(mcmcIte+2-burnin), "\n", file=outfile,append=T);
    if (header)
      cat("Temperature","model_idx", dataObj$xNames, "xnum", "freq", "\n", file=outfile,append=T);
    
    cumFreq <- cumsum(freqChain[1:currentPos]);
    startPos <- which(cumFreq>=burnin)[1];
    sFreq <- cumFreq[startPos]-burnin+1;  
    visited.cnt <- tapply(c(sFreq, freqChain[(startPos+1):currentPos]), gammaChain[startPos:currentPos], sum );
    visited.cnt <- sort(visited.cnt[1:min(length(visited.cnt),keepN )], decreasing=T);

    mapply(function(gm,fq) {gM <- gammaN2M(gm);
                            cat(temperature,gm, gM, sum(gM), fq,"\n", file=outfile ,append=T)},
           as.numeric(names(visited.cnt)),
           visited.cnt
           );
  } 
  
  h1 <- priorPara[1]
  tau11 <-  priorPara[2]
  if (length(priorPara)==2) {
    tau00 <- priorPara[2] 
  } else {
    tau00 <- priorPara[3] 
  }
  
  nNodes <- attr(edges,'nNodes')
  neighbors <- attr(edges,'neighbors')
  w <- attr(edges,'w')

  
  if (sampleFrmPrior) {
    onestepFun <- Gibbs.onestep.sampleFrmPrior
    
    currentNodesLabel <- if(is.null(initLabel) || is.na(initLabel[1]) ) 2*rbinom(nNodes,1,0.5)-1  else initLabel 
  }
  else {
    onestepFun <- Gibbs.onestep
    lik.1 <- likpw(nodesData, rep(1,nNodes), priorPara2)
    lik.0 <- likpw(nodesData, rep(-1,nNodes))
    likRatio <- (lik.1/lik.0)^(1/temperature) 
    
    currentNodesLabel <- if(is.null(initLabel)) 2*rbinom(nNodes,1,0.5)-1  else if(is.na(initLabel[1])) 2*(likRatio>1) - 1 else initLabel 
  }
  
  
  currentPos <- 1
  
  allNodesLabel <- matrix(0, ncol=nNodes, nrow=iteNum+1)  
  allPostLik <- rep(0,iteNum+1)
  
  
  allNodesLabel[1,] <- currentNodesLabel
  if (sampleFrmPrior)   { 
    allPostLik[1] <- graphPrior(currentNodesLabel,edges,h1, tau00, tau11) 
  }
  else {
    allPostLik[1] <- graphPrior(currentNodesLabel,edges,h1, tau00, tau11)*prod( ifelse( currentNodesLabel==1,lik.1,lik.0 ))
  }
  
  retVal <- list(env=environment(onestepFun), onestepFun=onestepFun,getChainLast=getChainLast,setChainLast=setChainLast)
  class(retVal) <- "Gibbs.t";
  return(retVal);
}


Gibbs.mstart <- function( mstart, iteNum, edges, priorPara, priorPara2=NULL,  initLabel=NULL, 
                          nodesData=NULL,sampleFrmPrior=F){
  temperature <- rep(1,mstart)
  allowSwap <- F
  
  if( !is.null(temperature)) {
    nT <- length(temperature);
  }
  else {
    temperature <- 1;
    nT <- 1
  }

  GibbsT <- vector("list",length=nT);
  
  for ( j in 1:nT ) {
    GibbsT[[j]] <- Gibbs.t(iteNum, edges, priorPara, priorPara2, temperature[j],initLabel=initLabel[[j]],
                           nodesData=nodesData, sampleFrmPrior=sampleFrmPrior);
  } 

  maxPost <- -Inf
  
  emcAcceptRate <- 0;
  mcmcPost <- vector(length=nT);
  for ( i in 1:iteNum){
    if ( i %%1000 ==0) print(i);
    for ( j in 1:nT ) {
       mcmcPost[j] <- GibbsT[[j]]$onestepFun();
       if(maxPost< mcmcPost[j]){
         maxPost <- mcmcPost[j]
         maxPostLab <- GibbsT[[j]]$getChainLast()
       }
     } 
    
    
    if (nT>1  && allowSwap) {
      for (j in 1:nT ) {
        
        xi <- sample(1:nT,1);
        xj <- if(xi==1) 2 else if(xi==nT) nT-1 else xi+sample(c(1,-1),1);
        probEx <- min(1, (mcmcPost[xi]/mcmcPost[xj])^(1/temperature[xj] - 1/temperature[xi]));
        if (rbinom(1, 1, probEx)) { 
          
          emcAcceptRate <- emcAcceptRate+1;
          xjLast <- GibbsT[[xj]]$getChainLast();
          GibbsT[[xj]]$setChainLast( GibbsT[[xi]]$getChainLast() );
          GibbsT[[xj]]$setChainLast( xjLast);          
        }
      } 
    } 
  } 
  
  if (nT>1  && allowSwap) {
    cat("Emc Accept Rate:", emcAcceptRate/nT/iteNum, 
      "\n");
  }
  attr(GibbsT,"maxPost") <- maxPost
  attr(GibbsT,"maxPostLab") <- maxPostLab
  return(GibbsT)
}



Gibbs.mcmc <- function( iteNum, edges, priorPara, priorPara2,  nodesData=NULL, initLabel=NULL,
                   sampleFrmPrior=F, returnPostLik=T, printI=T){ 
  if(!sampleFrmPrior && is.null(nodesData)) stop("Observed data for each node is required!")
  
  
  
  MRF.condProb <- function(nodei,nodesLabel) { 
    i.neighbor <- neighbors[[nodei]] 
    i.neighbor.0 <- i.neighbor[nodesLabel[i.neighbor]==-1]
    i.neighbor.1 <- i.neighbor[nodesLabel[i.neighbor]== 1]
    
    
    logitp <- h1 + tau11*( w[nodei]*sum(nodesLabel[i.neighbor]== 1)+sum(w[ i.neighbor.1])) -
                   tau00*( w[nodei]*sum(nodesLabel[i.neighbor]==-1)+sum(w[ i.neighbor.0]))   
    1/(1+exp(-logitp))
  }
  
  
  Gibbs.onestep.sampleFrmPrior <- function(nodesLabel) {
    nextNodesLabel <- nodesLabel
    for (nodei in 1:nNodes) {
      nextNodesLabel[nodei] <- 2*rbinom(1,1, MRF.condProb(nodei, nextNodesLabel))-1
    }
    nextNodesLabel
  }
  
  Gibbs.onestep <- function(nodesLabel) {
    nextNodesLabel <- nodesLabel
    for (nodei in 1:nNodes) {
      condProb <- MRF.condProb(nodei, nextNodesLabel)
      if(truncCondProb) condProb <- min(max(condProb,truncCondProbRange[1]),truncCondProbRange[2]) 

      prob <- (likRatio[nodei]*condProb) / (likRatio[nodei]*condProb + (1-condProb))
      nextNodesLabel[nodei] <- 2*rbinom(1,1, prob)-1
    }
    nextNodesLabel
  }
  h1 <- priorPara[1]
  tau11 <-  priorPara[2]
  if (length(priorPara)==2) {
    tau00 <- priorPara[2] 
  } else {
    tau00 <- priorPara[3] 
  }
  
  nNodes <- attr(edges,'nNodes')
  neighbors <- attr(edges,'neighbors')
  w <- attr(edges,'w')

  
  if (sampleFrmPrior) {
    Gibbs.onestep.fun <- Gibbs.onestep.sampleFrmPrior
    
    currentNodesLabel <- if(is.null(initLabel) || is.na(initLabel[1])) 2*rbinom(nNodes,1,0.5)-1  else initLabel 

  }
  else {
    Gibbs.onestep.fun <- Gibbs.onestep
    likRatio <- likpw(nodesData, rep(1,nNodes), priorPara2)/likpw(nodesData, rep(-1,nNodes)) 
    
    currentNodesLabel <- if(is.null(initLabel)) 2*rbinom(nNodes,1,0.5)-1  else if(is.na(initLabel[1])) 2*(likRatio>1) - 1 else initLabel
  }
  
  
  allNodesLabel <- matrix(0, ncol=nNodes, nrow=iteNum)  

  for (i in 1:iteNum) {
    if (i %% 5000==0 && printI) print( i );
    
    currentNodesLabel <- Gibbs.onestep.fun(currentNodesLabel)
    allNodesLabel[i,] <- currentNodesLabel
  }
  list(allNodesLabel=allNodesLabel)
}

icm <- function( nodesData, edges, priorPara, priorPara2,  initLabel=NULL, 
                searchMethod=c("random", "DFS","BFS"), printI=T, updatemethod=1) { 
   MRF.condProb <- function(nodei,nodesLabel) { 
    i.neighbor <- neighbors[[nodei]] 
    i.neighbor.0 <- i.neighbor[nodesLabel[i.neighbor]==-1]
    i.neighbor.1 <- i.neighbor[nodesLabel[i.neighbor]== 1]
    
    
    logitp <-  h1 + tau11*( w[nodei]*sum(nodesLabel[i.neighbor]== 1)+sum(w[ i.neighbor.1])) -
                    tau00*( w[nodei]*sum(nodesLabel[i.neighbor]==-1)+sum(w[ i.neighbor.0]))   
    1/(1+exp(-logitp))
  }
  
  updateLabel <- function(nodei) { 
    
    p1.nodei <- MRF.condProb(nodei,currentNodesLabel)
    currentNodesLabel[nodei] <<- 2*((likRatio[nodei] * (p1.nodei /(1-p1.nodei) ))>1 )-1
    
  }
  
  
  updateLabel.allnei <- function(nodei) { 
    
    p1.nodei <- MRF.condProb(nodei,currentNodesLabel)
    if(truncCondProb) p1.nodei <- min(max(p1.nodei,truncCondProbRange[1]),truncCondProbRange[2]) 

    testNodesLabel.0 <- testNodesLabel.1 <- currentNodesLabel
    testNodesLabel.1[nodei] <-  1
    testNodesLabel.0[nodei] <- -1
    pseudolik.ratio <- likRatio[nodei] * (p1.nodei /(1-p1.nodei))
    for ( j in neighbors[[nodei]] ){
      aa <- MRF.condProb(j,testNodesLabel.1)
      bb <- MRF.condProb(j,testNodesLabel.0)
      if(truncCondProb) {
        aa <- min(max(aa,truncCondProbRange[1]),truncCondProbRange[2]) 
        bb <- min(max(bb,truncCondProbRange[1]),truncCondProbRange[2]) 
      }

      pj.ratio <- if (currentNodesLabel[j]==1) aa/bb else (1-aa)/(1-bb)
      pseudolik.ratio <- pseudolik.ratio * pj.ratio
    }
    currentNodesLabel[nodei] <<- 2*(pseudolik.ratio >1)-1
    
    
    if(is.na(currentNodesLabel[nodei])) browser()
  }
  
  dfs <- function(initNode){ 
    visited[[initNode]] <<- T;
    for(child in neighbors[[initNode]]) {
      if(!visited[[child]] && child!=initNode) dfs(child);
    }
    
    update.fun(initNode);
  }
  
  bfs <- function(initNode) {
    visited[[initNode]] <<- T;
    
    update.fun(initNode);
    for(child in neighbors[[initNode]]) {
      if(!visited[[child]] && child!=initNode ) dfs(child);
    }
  }
  randsearch <- function(initNode) {
    for(nodeInd in sample(nNodes) ){
      update.fun(nodeInd)
    }
  }
  update.fun <- if(updatemethod==1) updateLabel else updateLabel.allnei
  h1 <- priorPara[1]
  tau11 <-  priorPara[2]
  if (length(priorPara)==2) {
    tau00 <- priorPara[2] 
  } else {
    tau00 <- priorPara[3] 
  }

  
  nNodes <- attr(edges,'nNodes')
  neighbors <- attr(edges,'neighbors')
  w <- attr(edges,'w')
  
  lik.1 <- likpw(nodesData, rep(1,nNodes), priorPara2)
  lik.0 <- likpw(nodesData, rep(-1,nNodes))
  likRatio <- lik.1/lik.0  
 
  currentNodesLabel <- if(is.null(initLabel)) 2*rbinom(nNodes,1,0.5)-1  else if(is.na(initLabel[1])) 2*(likRatio>1) - 1 else initLabel 
  searchFun <- if(searchMethod[1]=="BFS") bfs else if(searchMethod[1]=="DFS") dfs else randsearch
  converge <- F
  iteNum <- 1
  while (! converge ) {
    if (printI)  print(iteNum)
    prevNodesLabel <- currentNodesLabel
    
    visited <- rep(F,nNodes) 
    
    searchFun(sample(1:nNodes,1))
    converge <- all(prevNodesLabel == currentNodesLabel)
    iteNum <- iteNum+1
  }
  attr(currentNodesLabel,'posterior') <- graphPrior(currentNodesLabel,edges,h1, tau00, tau11)*
    prod( ifelse( currentNodesLabel==1,lik.1,lik.0 ))
  currentNodesLabel
}

logistic.estpara <- function(labelconfig,nodeNeighbors,w, tau.equal=T,lambda=0.5,maxite=100) { 
  if (tau.equal) {
    x <- cbind(x0=1,x1=mapply(function(neiIdx,wi) wi*sum(labelconfig[neiIdx])+ sum(w[neiIdx]*labelconfig[neiIdx]), nodeNeighbors, w))
    beta0 <- c(0,0) 
  } else {
    x <- cbind(x0=1,x1=mapply(function(neiIdx,wi)  sum(wi+w[neiIdx][labelconfig[neiIdx]== 1]),nodeNeighbors,w),
                    x2=mapply(function(neiIdx,wi) -sum(wi+w[neiIdx][labelconfig[neiIdx]==-1]),nodeNeighbors,w)) 
    beta0 <- c(0,0,0)
  }
  y <- (labelconfig+1)/2
  np <- ncol(x)
  converg <- rep(F,np)
  ite <- 1
  while(!all(converg) && ite<maxite) {
    ite <- ite+1
    linearpred <- as.numeric(x %*% beta0)
    prob <- 1/(1+exp(- linearpred))
    
    beta. <- beta0 + solve(t(x) %*% diag(prob*(1-prob)) %*% x + diag(2*lambda,np) ) %*% ( t(x)%*% (y-prob )-2*lambda*beta0 )
    
    converg <- abs(beta.-beta0)/(0.1+abs(beta.)) < 1e-4
    
    beta0 <- beta.
  }
  if(!all(converg)) {
    cat("Ridge logistic regression not converge after", maxite,"iterations! convergence of b_0...,b_p:", converg, "\n");
  }
  attr(beta.,'convergence') <- converg
  
  
  beta.
}


est.tau.h <- function( nodesData, edges, tau.equal=T, initPara=NULL, initLabel=NULL,lambdalist=0.5,maxite=200,
                      paraestonly=T) {
  
  nNodes <- attr(edges,'nNodes')
  neighbors <- attr(edges,'neighbors')
  w <- attr(edges,'w')

  
  estimateMuSiqsq <- function(nLabel) { 
    tempd <- y[nLabel==1]
    if((nn <- length(tempd))==0) c(2,1) else c(mean(tempd), if(nn==1) 1 else max(var(tempd)*(nn-1)/nn -1,0) ) 
  }
  
  if(tau.equal && length(initPara)==3) stop("Number of prior parameters should be 2!");
  nodeD.pval <- attr(nodesData,'p')

  
  nodesLabel <- if(is.null(initLabel)||is.na(initLabel[1])) (nodeD.pval<0.1)*2-1  else initLabel 

  if(class(nodesData)=="pval") {
    y <-  attr(nodesData,'y')
    
    priorPara2 <- estPriorPara2(nodesData)
  }

  retPara <- retPara2 <- crossValidLoss <- allNodesLabel <- NULL
  
  for (lambda in lambdalist) {
    priorPara <- if(is.null(initPara)) logistic.estpara(nodesLabel,neighbors,w, tau.equal=tau.equal,lambda=lambda) else initPara

    converg <- rep(F,length(priorPara))
    ite <- 0
    while(!all(converg) && ite<maxite ) {
      ite <- ite+1

      
      nodesLabel <- icm( nodesData, edges, priorPara, priorPara2, initLabel=nodesLabel, searchMethod=c("DFS","BFS","random"), printI=F,updatemethod=2 )
      prePriorPara <- priorPara
      priorPara <- logistic.estpara(nodesLabel,neighbors,w, tau.equal=tau.equal,lambda=lambda)

      converg <- abs(priorPara-prePriorPara)/(0.1+abs(priorPara))<1e-3
    }
    if(!all(converg)) cat("Lambda=",lambda,": ICM not converge after", maxite,"iterations!", "\n");
    
    retPara <-  cbind(retPara, as.numeric(priorPara ))
    retPara2 <- cbind(retPara2,priorPara2)
    allNodesLabel <- cbind( allNodesLabel, nodesLabel)
    
    
    crossValidLoss <- append(crossValidLoss, mean( abs( (nodeD.pval<=0.2) - (nodesLabel==1)) ))
  }
  rownames(retPara) <- c("h1", "tau11","tau00")[1:(nrow(retPara))]
  colnames(retPara) <- colnames(retPara2) <- lambdalist
  if(paraestonly) {
    list(paraest=retPara,paraest2=retPara2,crossValidLoss=crossValidLoss,allNodesLabel=allNodesLabel)}
  else {
    list(paraest=retPara,paraest2=retPara2,crossValidLoss=crossValidLoss,allNodesLabel=allNodesLabel,converge=converg)
  }
}



dispRst <- function(truenode,trueEdge,rst) {
  cat("True Node Value",truenode,sep=",","\n")

  netNeighbour <- rep(0, length(truenode))
  nu <- apply(trueEdge, 1, function(edge){
    if(truenode[edge[1]]==truenode[edge[2]]) deg <- 1 else deg <- -1;
    netNeighbour[edge[1]] <<- netNeighbour[edge[1]]+deg;
    netNeighbour[edge[2]] <<- netNeighbour[edge[2]]+deg;
  } )

  cat("Net Same Neighbourhood",netNeighbour,sep=",","\n" )
  cat("Reject Null by Z test,")
  cat(rst$reject, sep=",", "\n")
  cat("Reject Null by Max Posterior prob. ,")
  cat(rst$rejectMaxPostLik, sep=",","\n")
  cat("Reject Null by Most Freqent,")
  cat(rst$rejectMostFreq,sep=",", "\n")
  cat("Reject Null by (posterior odds)>1,")
  cat(apply(rst$postRatio>1,  2,mean),sep=",", "\n")
  cat("Reject Null by (posterior odds)>2,")
  
  cat(apply(rst$postRatio>2,  2,mean),sep=",", "\n")
  cat("Reject Null by (posterior odds)>2.5,")
  cat(apply(rst$postRatio>2.5,2,mean),sep=",", "\n")
  cat("Reject Null by (posterior odds)>3,")
  cat(apply(rst$postRatio>3,  2,mean),sep=",", "\n")
  cat("Reject Null by (posterior odds)>3.5,")
  cat(apply(rst$postRatio>3.5,2,mean),sep=",", "\n")
  cat("Reject Null by (posterior odds)>4,")
  cat(apply(rst$postRatio>4,  2,mean),sep=",", "\n")
  cat("Avg. posterior odds,")
  cat(rst$meanPostRatio,sep=",", "\n");
  cat("Avg. Poterior prob.,")
  cat(rst$meanAvgPost,sep=",", "\n")
}


simuTrueLab.Power <- function(nodesAsso, trueEdges, testRunNum, mcmcIte, mcmcBurnin, nNodes=length(nodesAsso),
                              classnm=c("pval","norm"), normal.std=1, priorPara=NULL,tau.equal=T, ridge.lambda=0.2) {
  nodesMean <- as.numeric(nodesAsso==1) 
  mcmcNoBurnin <- mcmcIte - mcmcBurnin
  
  
  reject <- rep(0,nNodes)
  priorPara.est <- matrix(0, ncol=5-tau.equal, nrow=testRunNum)
  maxPostLik <- mostFreq <- postRatio <- avgPost <- matrix(0, ncol=nNodes, nrow=testRunNum) 
  for (trn in 1:testRunNum) {
    if (trn%%100==0) cat("Running ", trn, "\n")
    
    nodesD <- simulateData(classnm[1], nodesMean, nodesStd=normal.std)
    
    reject <- reject + (attr(nodesD,"p")<0.05)

    
    prior.est <- if(is.null(priorPara)) est.tau.h(nodesD,trueEdges, tau.equal=tau.equal,lambda=ridge.lambda) else priorPara
    priorPara.est[trn,] <- prior.est

    
    rst <- Gibbs.mcmc( mcmcIte,  trueEdges, priorPara=prior.est, nodesData=nodesD, sampleFrmPrior=F, printI=F)
    
    
    maxPostLikInd <- which.max(rst$allPostLik) 
    maxPostLik[trn,] <- rst$allNodesLabel[maxPostLikInd,]

    allLabelNoBurnin <- rst$allNodesLabel[(mcmcBurnin+1):mcmcIte,]
    
    labelNum <- apply(allLabelNoBurnin, 1, label2num)
    mostFreq[trn,] <- num2label(as.numeric(names(which.max(table (labelNum)))),nNodes)

    
    cntOnes <- apply(allLabelNoBurnin==1, 2, sum)
    postRatio[trn,] <- cntOnes/(mcmcNoBurnin-cntOnes)

    
    avgPost[trn,] <- cntOnes/mcmcNoBurnin
  }
  colnames(priorPara.est) <- names(prior.est)
  retval <- list(reject=reject/testRunNum,
       rejectMaxPostLik=apply(maxPostLik==1,2,mean),
       rejectMostFreq=apply(mostFreq==1,2,mean),
       rejectPostRatio=apply(postRatio>2, 2, mean),  
       meanPostRatio=apply(postRatio,2,mean),
       meanAvgPost=apply(avgPost,2,mean),
       maxPostLik=maxPostLik,
       mostFreq=mostFreq,
       postRatio=postRatio,
       avgPost=avgPost,
       priorPara.est=priorPara.est)
  dispRst(nodesAsso,trueEdges,retval)
  retval
}


dispSamplePrior <- function(rst,nNodes,mcmcIte=nrow(rst$allNodesLabel),mcmcBurnin=0) {
  mcmcNoBurnin <- mcmcIte - mcmcBurnin
  priorLabelDraw <- rst$allNodesLabel[(mcmcBurnin+1):mcmcIte, ]

  
  maxPriorProbInd <- which.max(rst$allPostLik) 
  maxPriorProb <- rst$allNodesLabel[maxPriorProbInd,]
  
  
  labelNum <- apply(priorLabelDraw, 1, label2num)
  mostFreq <- num2label(as.numeric(names(which.max(table (labelNum)))),nNodes)

  
  cntOnes <- apply(priorLabelDraw==1, 2, sum)
  postRatio <- cntOnes/(mcmcNoBurnin-cntOnes)

  
  avgPriorProb <- cntOnes/mcmcNoBurnin

  cat("Max. prior probability,", paste(maxPriorProb, collapse=","), "\n");
  cat("Most freqently:,", paste(mostFreq, collapse=","), "\n");
  cat("Avg. prior probability,", paste(avgPriorProb, collapse=","), "\n")
  cat("Correlation of prior probability,\n")
  print( cor(priorLabelDraw))
  cat( "\n")

}

samplePrior <- function(trueEdges, nNodes, priorPara, mcmcIte, mcmcBurnin ){
  mcmcNoBurnin <- mcmcIte - mcmcBurnin
  
  rst <- Gibbs.mcmc( mcmcIte,  trueEdges, priorPara,  nNodes=nNodes, sampleFrmPrior=T, returnPostLik=T, printI=T)
  dispSamplePrior(rst,nNodes, mcmcIte, mcmcBurnin)
  rst
}



hist.nh <- function(cnt1,nNodes,cnt11,nEdges, priorPara=NULL,extrattl="") {
  if(is.null(priorPara)) {
    mainttl <- extrattl
  }
  else {
    priorPara <- signif(priorPara,3)
    h1 <- priorPara[1];tau11 <-  priorPara[2];tau00 <- priorPara[3] 
    mainttl <-  if(extrattl=="") bquote(list(h==.(h1), tau[1]==.(tau11), tau[0]==.(tau00))) else bquote(
                                                     list(.(extrattl), h==.(h1), tau[1]==.(tau11), tau[0]==.(tau00)))
  }
  hist(cnt1,xlim=c(0,nNodes), xlab="# nodes 1", main= mainttl)
  legend("topright",legend=bquote("Avg.="* .(round(mean(cnt1),2))/.(nNodes)==.(round(mean(cnt1)/nNodes,3))))
  
  hist(cnt11,xlim=c(0,nEdges), xlab="# edges 1-1", main= mainttl)
  legend("topright",legend=bquote("Avg.="* .(round(mean(cnt11),2))/.(nEdges)==.(round(mean(cnt11)/nEdges,3))))
}

estPriorPara.samplePrior <- function(trueEdges,  h1vec, tau00vec, tau11vec=NULL, mstart,mcmcIte, mcmcBurnin ){
  nNodes <- attr(trueEdges,'nNodes')
  nEdges <- nrow(trueEdges)
  
  aacc <- acnt1 <- acnt11 <- acnt00 <- apmean <- NULL
  mcmcNoBurnin <- mcmcIte - mcmcBurnin
  for(h1 in h1vec) {
    for(tau00 in tau00vec){
      alltau11 <- if(is.null(tau11vec)) tau00 else tau11vec
      for(tau11 in alltau11){
        priorPara <- c(h1,tau11,tau00)
        rst <- Gibbs.mstart( mstart, mcmcIte,  trueEdges, priorPara, sampleFrmPrior=T)
        ana <- gmstart.rst(trueEdges,rst,mcmcBurnin)
        apmean <- cbind(apmean, ana$postMean)
        aacc <- cbind(aacc,ana$acceptT)  
        acnt1 <-  cbind(acnt1, as.numeric(ana$cnt1))
        acnt11 <- cbind(acnt11,as.numeric(ana$cnt11))
        acnt00 <- cbind(acnt00,as.numeric(ana$cnt00))
      } 
    }
  }
  list(acc=aacc, pmean=apmean, cnt1=acnt1, cnt11=acnt11, cnt00=acnt00)
}



g1start.rst <- function(edges, rst1start,mcmcBurnin) { 
  mcmc1start <- get("allNodesLabel", env=rst1start$env)
  mcmcNum <- nrow(mcmc1start)
  
  
  cnt1 <- apply(mcmc1start==1, 1, sum)

  
  cnt11 <- apply(mcmc1start, 1, function(lab) sum(lab[edges[,1]] + lab[edges[,2]]==  2))
  cnt00 <- apply(mcmc1start, 1, function(lab) sum(lab[edges[,1]] + lab[edges[,2]]== -2))

  
  if (mcmcNum>mcmcBurnin) {
    mcmc1start <- mcmc1start[(mcmcBurnin+1):mcmcNum,]
    mcmcNum <- nrow(mcmc1start)
  }
  
  
  acceptT <- 0
  if (mcmcNum>=2) {
    for (m in 2:mcmcNum) {
      acceptT <- acceptT + all(mcmc1start[m,]==mcmc1start[m-1,])
    }
  }
  allLab02 <- apply(mcmc1start, 1,function(rowlab) paste(rowlab+1,collapse=""))
  visited <- as.data.frame(table(allLab02)) 
  vFreq <- nrow(visited)

  
  iMostFreq <- which.max(visited$Freq)
  mostFreq <- as.numeric(strsplit(as.character(visited[iMostFreq,"allLab02"]),"")[[1]])-1
  
  margSum <- apply( mcmc1start==1, 2, sum)

  return(list(mcmcNum=mcmcNum, acceptT=acceptT, cnt1=cnt1, cnt11=cnt11, cnt00=cnt00,visited=visited,
              
              vFreq=vFreq, mostFreq=mostFreq,
              margSum=margSum))
}

gmstart.rst <- function(edges, rstmstart,mcmcBurnin) { 
  map <- attr(rstmstart,"maxPostLab")$nodeLabel
  maxPost <- attr(rstmstart,"maxPost")
  
  postMargSum <- postFreqLab02 <- NULL
  mcmcNum <- acceptT <- rep(0,length(rstmstart))
  cnt1 <- cnt11 <- cnt00 <- NULL
  for (i in 1:length(rstmstart)) {
    if(is.null(rstmstart[[i]])) next

    rst1s <- g1start.rst(edges, rstmstart[[i]], mcmcBurnin)
    mcmcNum[i] <- rst1s$mcmcNum
    acceptT[i] <- rst1s$acceptT
    cnt1 <- cbind(cnt1,rst1s$cnt1)
    cnt11 <- cbind(cnt11,rst1s$cnt11)
    cnt00 <- cbind(cnt00, rst1s$cnt00)
    
    postMargSum <- rbind(postMargSum, rst1s$margSum)
    
    visited <- rst1s$visited
    colnames(visited)[2] <- paste("Freq",i,sep="")
    postFreqLab02 <- if(i==1) visited else merge(postFreqLab02,visited, all=T)
  }

  
  postMean <- apply( postMargSum, 2, sum)/sum(mcmcNum)
  
  totFreq <- apply(postFreqLab02[,-1,drop=F],1,sum,na.rm=T)
  indMostFreq <- which.max(totFreq)
  mostFreq <- as.numeric(strsplit(as.character(visited[indMostFreq,1]),"")[[1]])-1
  list(accept=sum(acceptT)/sum(mcmcNum-1), acceptT=acceptT/(mcmcNum-1), mcmcNum=sum(mcmcNum), postMean=postMean, mostFreq=mostFreq,
       mostFreqProp=totFreq[indMostFreq]/sum(mcmcNum), cnt1=cnt1, cnt11=cnt11,cnt00=cnt00)
}

