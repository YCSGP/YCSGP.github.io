## Data Pre-processing for Label-free Multiple Reaction Monitoring (MRM) Experiments
## Lisa M Chung, Christopher M Colangelo, and Hongyu Zhao
## Preparation for the submission of MDPI Biology, Special Issue: Advances in Proteomics Methods


## Required Data Format ###################################################
## 1) "transition.info" : Data Frame of Protein, Peptide, and Transition                              ##
##   - These three variables should be character strings format.                                        ##
##   - Number of rows should be the same as the number of transitions.                           ##
##                                                                                                                                        ##
##   - Example:                                                                                                                    ##       
##   > transition.info[1:5, ]                                                                                                    ##
##                   Protein                   Peptide   Transition                                                        ##
##   1 ADT1_MOUSE   DFLAGGIAAAVSK       y(10)                                                           ##
##   2 ADT1_MOUSE   DFLAGGIAAAVSK        y(9)                                                            ##
##   3 ADT1_MOUSE   DFLAGGIAAAVSK        y(6)                                                            ##
##   4 ADT1_MOUSE   DFLAGGIAAAVSK        y(8)                                                            ##
##   5 ADT1_MOUSE   DFLAGGIAAAVSK      y(11)                                                            ##
##                                                                                                                                        ##
##                                                                                                                                        ##
## 2) "Area" and "RT": A matrix of peak intensities and                                                      ##
##                                observed retention times from MRM experiments, respectively.     ##
##    - Each row should correspond to each transition in transition.info.                              ##
##    - Each column indicates each MRM run.                                                                     ##
##                                                                                                                                        ##
##    > Area[ 1:5 , 1:3 ]                                                                                                        ##
##                  1_1              1_2               1_3                                                                       ##      
## 1   3991436.2    3426231.2    3603866.5                                                                       ##
## 2   4010127.6    4069327.6    4325630.7                                                                       ##
## 3   3013942.5    3011190.8    2824829.7                                                                       ##
## 4   1245235.9    1269667.3    1233300.9                                                                       ##
## 5     654219.3      641204.1      635197.4                                                                       ##
########################################################################


## Required Packages ######################################################
library(lattice)
library(MASS)
library(affy)
library(limma)

########################################################################
## I. Data Quality Assessment ################################################
########################################################################

## 1) Observed retention time deviation and coefficient of variation ###################

## transition-specific median retention time
med.RT <- apply(RT, 1, median, na.rm = TRUE) 

## retention time deviation from median
dev.RT <- RT - med.RT

## construct data frame for graphic
df.RT <- stack(as.data.frame(dev.RT))
df.RT$med.RT <- rep(med.RT, ncol = dev.RT)

## graphic (RT deviation graph as in  Figure 1 (a))
cols <- rainbow(ncol(RT)+3)[1:ncol(RT)]
xyplot( values ~ med.RT, groups = ind,  col = cols, df.RT, 
		   type = c("g", "p"), xlab = "Median RT (mins) for each fragment",
		   ylab = "RT deviation from the median",
		   auto.key = list(space = "right", col = cols, points = list(pch = "")))
   

## Adjusted Coefficient of Variation
med.factor <- apply(RT, 2, median, na.rm = TRUE) -  median(RT, na.rm =TRUE)
adj.RT <- t( t(RT)  - med.factor )

adj.RT.CV <- apply(adj.RT, 1, function(x){sd(x, na.rm = TRUE)/mean(x, na.rm = TRUE)})


## 2) Weight calculation to assess the quality of each peak  ##################
##    ( Next four lines are to convert the data format. )
log2.area <- stack(as.data.frame( log2(Area) ) )
log2.area$Sample <- as.character(log2.area$ind)
log2.area$Protein <- rep(transition.info$Protein, ncol(Area))
log2.area$Peptide <- rep(transition.info$Peptide, ncol(Area))
log2.area$Transition <- rep(transition.info$Transition, ncol(Area))

log2.area$wts <- rep(NA, nrow(log2.area))

pep.id <- unique(log2.area$Peptide)
n.peps <- length(pep.id)

for(i in (1:n.peps)){
	
	sub <- subset( log2.area, Peptide == pep.id[i] )
	sub <- subset(sub, !is.na(values))
	
	
	m1 <- rlm( values ~ as.factor(Sample) + as.factor(Transition), sub, maxit = 50)
	
	sub$wts <- m1$w
	log2.area[ rownames(sub) , ]$wts <- sub$wts
		
}

## This is a weight matrix for each transition-level peak 
wts.matrix <- matrix(log2.area$wts, nr = nrow(Area), nc = ncol(Area),
								dimnames = list(rownames(Area), paste(colnames(Area), "weights", sep= "_")))

## 3) summarized for each transition ###########################
med.wts <- apply(wts.matrix, 1, median, na.rm = TRUE)
q1.wts <- apply(wts.matrix, 1, quantile, probs = 0.25, na.rm = TRUE)
q3.wts <- apply(wts.matrix, 1, quantile, probs = 0.75, na.rm  = TRUE) 
iqr.wts <- q3.wts - q1.wts

## Result Table 
QA.result <- data.frame( transition.info, 
										"adjusted.RT.CV"  = adj.RT.CV,
										"Median.Fragment.Weights" = med.wts,
										"IQR.Fragment.Weights" = iqr.wts,
										wts.matrix)

write.table(QA.result, "RT_CV_and_Quality_Weights.txt", 
					sep = "\t", col.names = TRUE, row.names= FALSE, quote= FALSE)

## 4) NUSE #############################################
## calculate summed weight for each peptide
WW <- apply(wts.matrix, 2, function(x, ind){ tapply(x, ind, sum, na.rm = TRUE) }, 
						ind = transition.info$Peptide)
	
WW <- 1/sqrt(WW) # total peptide weight

med.WW <- apply(WW2, 1, median, na.rm = TRUE)

nuse <- (WW)/med.WW ## NUSE

NUSE.summary <- data.frame( "Median" = apply(nuse, 2, median),
												 "IQR" = apply(nuse, 2, 
												 function(x){ quantile(x, 0.75) - quantile(x, 0.25)}))
												 
												 
write.table(NUSE.summary, "NUSE_for_sample_quality.txt", 
					sep = "\t", col.names = TRUE, row.names= TRUE, quote= FALSE)

#########################################################
##  II. Data Normalization ####################################
#########################################################
log2.area.matrix <- log2(Area) 

# if one would like to use transitions with adjusted CV < 5%
# log2.area.matrix <- log2.area.matrix[ adj.RT.CV < 0.05, ]

# if one would line to use transitions with weight median  > 0.9
# log2.area.matrix <- log2.area.matrix[ med.wts > 0.9, ]


## 1) Global median adjustment
med.factor <- apply(log2.area.matrix, 2, median, na.rm = TRUE)  - median(log2.area.matrix, na.rm = TRUE)
normed.area <- t( t(log2.area.matrix) - med.factor )

## 2) Quantile normalization
normed.area <- normalize.quantiles( log2.area.matrix )

## 3) Cyclic loess
normed.area <- normalizeCyclicLoess( log2.area.matrix )

## 4) Subset-based median adjustment (example)
subset.ind <- which(transition.info$Protein == "ADT1_MOUSE")
med.factor <- apply(log2.area.matrix[subset.ind, ], 2, median, na.rm = TRUE)  - 
						median(log2.area.matrix[subset.ind, ], na.rm = TRUE)
	
normed.area <-  t( t(log2.area.matrix) - med.factor )			
						
## 5) Invariant-set normalization
normed.area <- log2.area.matrix
normed.area[,] <- NA
area.ref <- apply(log2.area.matrix, 1, median, na.rm = TRUE)
inv.set <- list()

for(j in 1:ncol(log2.area.matrix)){

	y1 <- log2.area.matrix[,j]
	ind.na <- is.na(y1)

	ref1 <- area.ref[ !ind.na ]
	yy1 <- y1[ !ind.na ]

	foo <- normalize.invariantset(data = yy1, ref = ref1, prd.td = c(0.003, 0.007))
	foo2 <- approx(foo$n.curve$y, foo$n.curve$x, xout = yy1, rule = 2)$y 

	inv.set[[j]] <- foo$i.set

	normed.area[ names(yy1), j ] <- foo2

}

normed.area.df <- data.frame(transition.info, normed.area)
write.table(normed.area.df, "Normalized_PeakArea.txt", 
					col.names = TRUE, row.names = FALSE, quote = FALSE, sep= "\t")
					
					

