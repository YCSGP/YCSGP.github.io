## Lisa Chung, Hongyu Zhao, et al. 
## Bayesian modeling for paired RNA-seq experiment
## This is a testing version..

## 2-gp log-normal model

# Files: pairedBayesRsource.R is the R script, source file and
#        pairedBayesExample2.Rdat is an example, simulated data set.

### DETAILS ON pairedBayesR.R ###############################################
# Required package: MCMCpach, pscl
# input: 
#   pre.yy : pre-treatment count data
#   post.yy : post-treatment count data
#   (columns of pre.yy and post.yy should match (from the same individual)
#   (rows of pre.yy and post.yy should match (from the same gene/transcript)
#   pre.NN : library sizes for pre-treatment data
#   post.NN: library sizes for post-treatment data
#   num.iter: number of total MCMC iterations
#   n.burnin: number of a burnin iterations
#   sigma.alpha and sigma.beta: standard deviation to 
#                               Metropolis-Hastings Samplling for
#                               baseline expression
#############################################################################


# Example:
# 
# load("pairedBayesExample2.Rdat")
# source("pairedBayesRsource.R")
#
# res <- runMCMC(pre.yy = pre.yy, post.yy=post.yy, pre.NN = pre.NN, post.NN = post.NN)
# res$postprob : a matrix of posterior probabilities (EE and DE)
# res$foldChange: a vector of fold change estimated from 
#                 this Baysian approach, logged scale 
