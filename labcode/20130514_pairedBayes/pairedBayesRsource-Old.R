## 2-gp log-normal model


require(MCMCpack)
require(pscl)



######################################################################
## Distributions to generate proposal (jump functions)              ##
######################################################################

runMCMC <- function(pre.yy, post.yy, pre.NN, post.NN, num.iter = 6000, n.burnin = 3000, sigma.alpha = 0.1, sigma.beta = 0.5){

	#max.cnt <- (apply(cbind(pre.yy, post.yy), 1, max))
	#max.ind <- which(max.cnt < quantile(max.cnt, prob = 0.95))

	#pre.NN <- round(colSums(pre.yy[max.ind, ])*(sum(pre.yy[,1])/sum(pre.yy[max.ind,1])))
	#post.NN <- round(colSums(post.yy[max.ind, ])*(sum(pre.yy[,1])/sum(pre.yy[max.ind,1])))

	


	sigma2.alpha <- sigma.alpha^2
	sigma2.beta <- sigma.beta^2


	## try MCMC


	if(any(c(ncol(pre.yy)!=ncol(post.yy), nrow(pre.yy)!=nrow(post.yy), length(pre.NN) != length(post.NN)))){stop("should be paired counts")}


	num.sample <- nrow(pre.yy)




	# each column of the pre.yy and post.yy: from the same individual

	# parameters to estimate
	# : alpha, beta from lambda ~ gamma(shape = alpha, rate = beta)
	# : alpha.pos, beta.pos from chi ~ gamma(shape = alpha.pos, rate = beta.pos) when overexpressed in the post trt
	# : alpha.neg, beta.neg from chi ~ gamma(shape = alpha.neg, rate = beta.neg) when overexpressed in the pre trt
	# assume flat priors for parameters..

	######################################################################
	## some functions for the calculation of acceptance probabilities..  #
	######################################################################
	alpha.fn <- function(alpha, pre.yy, post.yy, 
		             pre.NN, post.NN, beta, chi){

		GG <- nrow(pre.yy); nn <- ncol(pre.yy)

		sum(lgamma(pre.yy + post.yy + alpha)) - 
        	alpha*sum(log(beta + outer(rep(1, GG), pre.NN) +  
		outer(chi, post.NN))) - 
		GG*nn*(lgamma(alpha) - alpha*log(beta)) 

	}


	alpha.g.fn <- function(alpha, pre.yy, post.yy, 
	             pre.NN, post.NN, beta, chi){

		# Here, alpha is a vector with length = G

		GG <- nrow(pre.yy); nn <- ncol(pre.yy)

	
		rowSums(lgamma(pre.yy + post.yy + outer(alpha, rep(1, nn) ))) - 
		nn*(lgamma(alpha) - alpha*log(beta)) - 		
        	alpha*rowSums(log(outer(beta, rep(1, nn)) + outer(rep(1, GG), pre.NN) + outer(chi, post.NN)))
	}




	beta.fn <- function(beta, pre.yy, post.yy, 
			    pre.NN, post.NN, alpha, chi){
	
		GG <- nrow(pre.yy); nn <- ncol(pre.yy)

		GG*nn*alpha*log(beta) -
			sum((pre.yy + post.yy + alpha)*
			log(beta + outer(rep(1, GG), pre.NN) +  
			outer(chi, post.NN)))
	}



	beta.g.fn <- function(beta, pre.yy, post.yy, 
			    pre.NN, post.NN, alpha, chi){

		# here, beta is a vector with length = G

		GG <- nrow(pre.yy); nn <- ncol(pre.yy)

		nn*alpha*log(beta) -
			rowSums((pre.yy + post.yy + alpha)*
			log(outer(beta, rep(1, nn)) + outer(rep(1, GG), pre.NN) +  
			outer(chi, post.NN)))
	}





	######################################################################
	## Define vectors/matrices for MCMC update                          ##
	######################################################################

	var.de.update <- mu.de.update <- var.ee.update <- rep(1, num.iter+1)

	chi.update <- zz.update <- alpha.update <- beta.update <- matrix(1, nr = nrow(pre.yy), nc = num.iter+1)
	pi.update <- matrix(0.5, nc = num.iter+1, nr = 2)



	######################################################################
	## Initialize..                                                     ##
	######################################################################


	pre.yy2 <- pre.yy; post.yy2 <- post.yy
	pre.yy2[pre.yy == 0] <- 0.5; post.yy2[post.yy == 0] <- 0.5


	d1 <- t(t(pre.yy2)/pre.NN); d2 <-  t(t(post.yy2)/post.NN); 
	diff <- 2^(rowMeans(log2(d2/d1)))

	diff[is.na(diff)] <- 1

	st <- rep(0, nrow(pre.yy))
	#st[800:1000] <- 1

	diff <- 0


	st <- sample(0:1, nrow(pre.yy), replace = TRUE)


	
	

	


	wts <- table(st)/length(st)

	#wts <- c(0.5, 0.5)

	# compute pdf from jump function: if z_g_pros (= u.zz) = +/- 1, sample chi_g from normal distribution..
	foo2 <- function(x,y, x.NN, y.NN){
	
		x <- t(t(x)/x.NN)
		y <- t(t(y)/y.NN)
			
		res <- log(x/y)
			
		res.mean <- apply(res, 1, FUN = function(x){median(x[!is.nan(x) & is.finite(x)])})
		res.var <-  apply(res, 1, FUN = function(x){var(x[!is.nan(x) & is.finite(x)])})

		res.rate <- res.mean/res.var
		res.shape <- res.rate * res.mean

		list("shape" = res.shape, "rate" = res.rate, "mean" = res.mean, "var" = res.var)
	}


	tmp <- foo2(x = post.yy2, y = pre.yy2, x.NN = post.NN, y.NN = pre.NN)



	proposal.mean <- tmp$mean
	proposal.var <- tmp$var

	proposal.var[proposal.var == 0] <- min(proposal.var[proposal.var > 0])
	#proposal.var[proposal.var > 10] <- 10
	#proposal.var[proposal.var > .5] <- .5

	#proposal.var[proposal.var > quantile(proposal.var, 0.8)] <- quantile(proposal.var, 0.8)
	#proposal.var[proposal.var < quantile(proposal.var, 0.4)] <- quantile(proposal.var, 0.4)


	proposal.var <- quantile(proposal.var, prob = 0.8)
	

	mm <- apply(pre.yy2, 1, function(x,n){mean(x/n)}, n = pre.NN)
	vv <- apply(pre.yy2, 1, function(x,n){var(x/n)}, n = pre.NN)

	mm[mm == 0] <- 1e-20
	vv[vv==0] <- 1e-20



	beta.update[,1] <- mm/vv
	alpha.update[,1] <- mm^2 / vv




	chi.update[,1] <- diff
	zz.update[,1] <- st
	pi.update[,1] <- wts


	## initialize hyper-paramters...
	var.de.update[1] <- 1

	mu.de.update[1] <- 0
	var.ee.update[1] <- 0.5

	######################################################################
	## Start the MCMC chain                                             ##
	######################################################################


	for(iter in 1:num.iter){
		# update alpha (Metropolis-Hastings)
		# proposal

		GG <- nrow(pre.yy)

		alpha.current <- alpha.update[,iter]
		beta.current <- beta.update[,iter]
		exp.chi.current <- exp(chi.update[,iter])

		## proposal: normal may not be reasonable, try beta or gamma...
		alpha.pros <- rnorm(GG, alpha.current, sigma.alpha) # normal
       

		u <- runif(GG)
		alpha.trans <- alpha.g.fn(alpha.pros, pre.yy, post.yy, 
         	                  pre.NN, post.NN, beta.current, exp.chi.current)-
          	               alpha.g.fn(alpha.current, pre.yy, post.yy, 
                 	          pre.NN, post.NN, beta.current, exp.chi.current)
		alpha.trans <- pmin(alpha.trans, 0)

		ind.alpha <- (log(u) < alpha.trans) & (alpha.pros > 1e-300)

		alpha.update[, iter+1] <- alpha.current
		alpha.update[ind.alpha, iter+1] <- alpha.pros[ind.alpha]
	

		# update beta  (Metropolis-Hastings)
		alpha.current <- alpha.update[iter+1]
		beta.pros <- rnorm(GG, beta.current, sigma.beta)
		#if(beta.pros <0){beta.pros > 1e-300}


		u <- runif(GG)

		beta.trans <- beta.g.fn(beta.pros, pre.yy, post.yy, 
                          pre.NN, post.NN, alpha.current, exp.chi.current)-
                         beta.g.fn(beta.current, pre.yy, post.yy, 
                          pre.NN, post.NN, alpha.current, exp.chi.current)
		beta.trans <- pmin(beta.trans, 0)

		ind.beta <- (log(u) < beta.trans) & (beta.pros > 1e-300)

		beta.update[, iter+1] <- beta.current
		beta.update[ind.beta, iter+1] <- beta.pros[ind.beta]



		# Jointly update (z_g and chi_g) 
        	# using Generalized Metropolis-Hastings
		beta.current <- beta.update[,iter+1]
		alpha.current <- alpha.update[,iter+1]



		pi.current <- pi.update[,iter]
		chi.current <- chi.update[,iter]
		zz.current <- zz.update[,iter]


		var.de.current <- var.de.update[iter]
		var.ee.current <- var.ee.update[iter]

		mu.de.current <- mu.de.update[iter]


  		# generate the z_g proposal
                u.zz <- runif(nrow(pre.yy))

                zz.pros <- sample(0:1, nrow(pre.yy), replace= TRUE)

                #zz.pros <- rep(0, length(u.zz))

                #zz.pros[which(u.zz > pi.current[1])] <- 1

                rr <- u.zz


                # generate chi_g given the z_g proposal

                proposal.mean.current <- (proposal.mean)
                proposal.mean.current[zz.current == 0] <- 0
                proposal.var.current <- proposal.var

                proposal.mean.pros <- proposal.mean
                proposal.mean.pros[zz.pros == 0] <- 0
                proposal.var.pros <- proposal.var


                chi.pros <- rnorm(nrow(pre.yy), mean = (proposal.mean.pros), sd =(sqrt(proposal.var.pros)))



   		## some quantities for acceptance probabilities
                gam.up.pros <- dnorm(chi.pros, mean = mu.de.current, sd = (sqrt(var.de.current)), log = TRUE)
                gam.up.current <- dnorm(chi.current, mean = mu.de.current, sd= (sqrt(var.de.current)), log = TRUE)

                gam.eq.pros <- dnorm(chi.pros, mean = 0, sd = (sqrt(var.ee.current)), log = TRUE)
                gam.eq.current <- dnorm(chi.current, mean = 0, sd = (sqrt(var.ee.current)), log = TRUE)




                norm.pros <- dnorm(chi.pros, mean = (proposal.mean.pros), sd = (sqrt(proposal.var.pros)), log = TRUE)
                norm.current <- dnorm(chi.current, mean = (proposal.mean.current), sd = (sqrt(proposal.var.current)), log = TRUE)




                log.trans.ratio <- (norm.pros) - (norm.current)


                log.fac.current <- post.yy * outer((chi.current), rep(1,ncol(post.yy))) - (pre.yy+post.yy+alpha.current) *log(beta.current + outer(rep(1, nrow(post.yy)), pre.NN) + outer(exp(chi.current), post.NN))



                log.fac.pros <- post.yy * outer((chi.pros), rep(1, ncol(post.yy))) - (pre.yy+post.yy+alpha.current) * log(beta.current + outer(rep(1, nrow(post.yy)), pre.NN) + outer(exp(chi.pros), post.NN))





                log.fac <- rowSums(log.fac.pros) - rowSums(log.fac.current)


                ind <- zz.current == 1 & zz.pros == 1
                rr[ind] <- (gam.up.pros[ind]) - (gam.up.current[ind])

                ind <- zz.current == 1 & zz.pros == 0
                rr[ind] <- (gam.eq.pros[ind]) - (gam.up.current[ind]) + log(pi.current[1]) - log(pi.current[2])


                ind <- zz.current == 0 & zz.pros == 1
                rr[ind] <- (gam.up.pros[ind]) - (gam.eq.current[ind])  + log(pi.current[2]) - log(pi.current[1])


                ind <- zz.current == 0 & zz.pros == 0
                rr[ind] <- (gam.eq.pros[ind]) - (gam.eq.current[ind])

                #rr2 <- exp(rr - log.trans.ratio)

                rr2 <- exp(log.fac + rr - log.trans.ratio )





		accept.prob <- pmin(rr2, 1)	
		does.accept <- runif(length(accept.prob)) < (accept.prob)
	
		zz.new <- zz.current
		zz.new[does.accept] <- zz.pros[does.accept]

		chi.new <- chi.current
		chi.new[does.accept] <- chi.pros[does.accept]

		chi.update[,iter+1] <- chi.new
		zz.update[,iter+1] <- zz.new


		# update weights (pi_, pi0, pi+)
		pi.update[,iter+1] <- rdirichlet(1, c(sum(zz.update[,iter+1] == 0)+1, sum(zz.update[,iter+1] == 1)+1))





		# update mu.de
		mu.de.pros <- rnorm(1, mean = sum((chi.new[zz.new == 1]))/sum(zz.new == 1), sd = sqrt(var.de.current/sum(zz.new == 1)))	

		mu.de.update[iter+1] <- mu.de.pros

		var.de.update[iter+1] <- rigamma(1, alpha = (sum(zz.new == 1)/2) , 
					beta = 0.5* sum( ((chi.new[zz.new == 1]) - mu.de.pros)^2 ))
	
		var.ee.update[iter+1] <- rigamma(1, alpha = (sum(zz.new == 0)/2) , 
					beta = 0.5 *  sum(((chi.new[zz.new == 0]))^2))

		if(iter %% 500 == 0){print(table(zz.new)); print(var.ee.update[iter+1]);  print(var.de.update[iter+1]);  print(mu.de.update[iter+1]); }





		if(iter %% 100 == 0){print(paste(iter, "iterations", Sys.time()))}

	}




	expr.class <- t(apply(zz.update[, n.burnin:num.iter], 1, 
			function(xx){c(mean(xx == 0), mean(xx == 1))}))

	colnames(expr.class) <- c("EE", "DE")

	chi.est <- rowMeans(chi.update[,n.burnin:num.iter])

	post2 <- apply(chi.update[,n.burnin:num.iter], 1, function(xx){mean(abs((xx)) > log(3))})


	list("postprob" = expr.class, "post2" = post2, "foldChange" = chi.est, "var.de" = var.de.update, "var.ee" =  var.ee.update, "mu.de" = mu.de.update, "mix.pi" = pi.update)

}