
analyze.iTRAQ <- function(){


	library(limma)

	normalize <- readline("Do you want to Normalize data (Y if Yes): \n")


	cat("Choose the data file: \n")

	file.nms <- file.choose()
	
	file.nms2 <- strsplit(basename(file.nms), ".txt")[[1]][1]
	file.nms2 <- strsplit(basename(file.nms2), ".tsv")[[1]][1]


	output.dir <- paste(  dirname(file.nms)  , "/iTRAQ_", Sys.Date(), "_", file.nms2, "/", sep = "" )
	dir.create( output.dir )

		
	fileout.nms <- file.nms2

	data.all <- read.delim(file.nms, stringsAsFactors = FALSE, na.strings = c("NA", "N/A"))


	cat("Only pick Annotation = auto and Cleavages = empty.. \n")
	data <- subset(data.all, Annotation == "auto ")
	data <- subset(data, Cleavages == "" | is.na(Cleavages))

	# check if it's 4-plex or 8-plex
	plex.index <- sapply( strsplit(colnames(data), "Area."), function(x){x[2]})
	plex.index <- plex.index[!is.na(plex.index)]

	plex <- length(plex.index)
	col.nms <- paste("Area.", plex.index, sep = "")

	if(plex != 4 & plex != 8){

		stop("The file is not from 4-plex nor 8-plex, or wrong column name. Please check.. \n")

	}


	area <- as.matrix(data[,col.nms])
	check.zero.area <- area == 0
	area[check.zero.area] <- NA

	log2.area <- log2(area)


	if (normalize == "Y" | normalize == "y"){

		cat("Normalizing.. \n")
		new.log2.area <- normalizeCyclicLoess(log2.area) 

		fileout.nms <- paste(output.dir, "AfterNormalization_", fileout.nms, ".tsv", sep = "")

	}else{

		cat("Just log2 transformation..")
		new.log2.area <- log2.area

		fileout.nms <- paste(output.dir, "NoNormalization_", fileout.nms, ".tsv", sep = "")


	}


	# Write a table (tsv format)
	results <- cbind(data[,is.na(sapply( strsplit(colnames(data), "Area."), function(x){x[2]}))], new.log2.area)
	write.table (results, fileout.nms, row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")


	# Now, let's make a graph..
	labs <- outer(plex.index, plex.index, paste)
	labs <- labs[upper.tri(labs)]


	# make a plot for raw data
	data.min <- min(as.vector(new.log2.area), na.rm = TRUE)
	data.max <- max(as.vector(new.log2.area), na.rm = TRUE)
	
	#x.range <- c(data.min, data.max )
	#y.range <- c( -data.max, data.max )


	cat("Making scatter plot for before/after normalization \n")
	pdf( paste(output.dir, "pairwiseScatterPlotForNormalization.pdf", sep = ""), height = plex*2, width = 8) 

	par(mfrow = c(plex/2, 2))


	for(ii in 1:length(labs)){
		
		new.ind <- which(colnames(new.log2.area) %in%  paste("Area.", strsplit(labs[ ii ], " ")[[1]], sep = ""))
		do.graph(log2.area, ind.x = new.ind[1], ind.y = new.ind[2], is.normalized = FALSE) 

		if(normalize == "Y" | normalize == "y"){

			do.graph(new.log2.area, ind.x = new.ind[1], ind.y = new.ind[2], is.normalized = TRUE) 

		}
	}

	dev.off()





	###################################################################################

	

	z <- readline("What cutoff do you want? \n")

	if(is.na(z) | is.null(z) | z == "" | is.na(as.numeric(z))){ z <- "1.5" }

	z <- log2(as.numeric(z))
	unlog.z <- 2^z
	
	
	choose.pair <- readline("Want to choose two groups to be compared? Not doing it will make ALL pairwise comparison. (y/n)")



	if(choose.pair == "y" | choose.pair == "Y"){

		
		comp1 <- readline("write the label names in GROUP1 to be compared, if more than one, separate by COMMAs: \n")
		comp2 <- readline("write the label names in GROUP2 to be compared, if more than one, separate by COMMAs: \n")



		comp1 <- paste("Area.", strsplit(comp1, ",")[[1]], sep = "")
		comp2 <- paste("Area.", strsplit(comp2, ",")[[1]], sep = "")

		protein.ids <- paste(results$Accession, results$Names, sep = "___")

		output.nms <- paste(output.dir, paste(comp1, collapse = ""), "vs", paste(comp2, collapse = ""), "_", unlog.z, "foldChange", Sys.Date(), sep = "")

		res <- protein.fc( xx = results[,comp1], yy = results[,comp2], protein.id = protein.ids, 
			cutoff = z, nms.xx = paste(comp1, ", "), nms.yy = paste(comp2, ", "), output.name = output.nms)  # 113,114 vs. 115,116


	}else{

		# This will make all pairwise comparison.
		for(i in 1:length(labs)){

 
			comp <- strsplit(labs[i], " ")[[1]]

			comp1 <- paste("Area.", strsplit(comp[1], ",")[[1]], sep = "")
			comp2 <- paste("Area.", strsplit(comp[2], ",")[[1]], sep = "")
			cat(paste(comp1 , "VS", comp2, ": \n"))


			protein.ids <- paste(results$Accession, results$Names, sep = "___")


			output.nms <- paste(output.dir, comp1, "vs", comp2, "_", unlog.z, "foldChange", Sys.Date(), sep = "")


			print(output.nms)

			res <- protein.fc( xx = results[,comp1], yy = results[,comp2], protein.id = protein.ids, 
							cutoff = z, output.name = output.nms, nms.xx = comp1, nms.yy = comp2)  

		}# end of for-loop for all pairwise comparison 

	}# end of doing all pairwise comparison
	







}










do.graph <- function(area, ind.x, ind.y, is.normalized, y.range, x.range){

	# This is a function for making a graph (MA plot like..)

	yy <- area[, ind.y]	
	xx <- area[, ind.x]

	ylab1 <- paste( colnames(area)[ind.y] , "-", colnames(area)[ind.x] )
	xlab1 <- paste( "avg between ", colnames(area)[ind.x] , "and" , colnames(area)[ind.y])
	main.title <- if(is.normalized){"Normalized"}else{"Raw data"}

		
	mm <- yy-xx
 	aa <- (xx+yy)/2
	plot(aa, mm, main = main.title, ylab = ylab1, xlab = xlab1, cex.lab = 1.2 )
 	abline(h=0, col = "grey")
 
	ll <- loess(mm ~ aa)
	points(ll$x, ll$fitted, col = "red", pch = ".", cex = 2)

}# end of do.graph function




###############################################################

# id <- "sp|Q3SZR3|A1AG_BOVIN___Alpha-1-acid glycoprotein OS=Bos taurus GN=ORM1 PE=2 SV=1" pval  0.5398


protein.fc <- function(xx, yy, protein.id, cutoff = NULL, output.name = NULL, nms.xx, nms.yy){ 

	xx <- as.matrix(xx)
	yy <- as.matrix(yy)

	xx2 <- rowMeans(xx)
	yy2 <- rowMeans(yy)



	fc <- rowMeans(yy) - rowMeans(xx)
	

	mm.fc <- tapply(fc, INDEX = protein.id, FUN = mean)
	sd.fc <- tapply(fc, INDEX = protein.id, FUN = sd)

	cv.fc <- sd.fc / mm.fc
	pval.fc <- tapply(fc, INDEX = protein.id, FUN = function(x){pval <- NA; if(length(x) >1 & sum(!is.na(x)) > 1){ pval <- t.test(x)$p.value}; pval} )


	num.peptides <- table(protein.id)

	median.fc <- tapply(fc, INDEX = protein.id, FUN = median)

	res <- cbind("unlog.fc" = 2^mm.fc, "avg.log2.fc" = mm.fc, "median.log2.fc" = median.fc, "sd.log2.fc" = sd.fc, "cv.log2.fc" = cv.fc, "p.val.from.ttest" = pval.fc, num.peptides)

	pick1.5 <- which( (abs(res[,3]) > cutoff & res[,7] > 2 & res[,6] < 0.1) |  (abs(res[,3]) > cutoff & res[,7] == 2 ))
	pick1.5.over2 <- which( (abs(res[,3]) > cutoff & res[,7] > 2 & res[,6] < 0.1) )
	

	



	#print(pick1.5)

	out <- res[pick1.5, ]

	if(length(pick1.5) == 1){ out <- as.matrix(t(out))}



	out.over2 <- out[ out[,7]  > 2, ] # pick a subset with proteins with >2 peptides

	if(sum(   out[,7]  > 2  ) == 1){ out.over2 <- as.matrix(t(out.over2))}



	

	protein.nms.picked <- rownames(res)[pick1.5]
	protein.nms.picked.over2 <- rownames(res)[pick1.5.over2]


	# do graphic and table if any over cutoff
	if(!is.null(out.over2) & nrow(out.over2) > 0 ){

		# Define color labels
		cols <- rep("red", nrow(out.over2))
		cols[(out.over2[,1] < 1)]<- "blue"

		pdf(paste(output.name, ".pdf", sep = ""), height = 5, width = 6)
		plot(xx2, yy2, main = paste("log2-scaled peptide-level plot \n on all proteins with", 2^cutoff, "fold change"), xlab = nms.xx,   ylab = nms.yy)
		abline(0,1, col = "grey")
		abline(cutoff, 1, col = "yellow")
		abline(-cutoff, 1, col = "yellow")


		for(kk in 1:length(protein.nms.picked.over2)){
			points(xx2[protein.id == protein.nms.picked.over2[kk]], yy2[protein.id == protein.nms.picked.over2[kk]], pch = 20, col = cols[kk])	
		}


		# Then, make individual ones
		for(kk in 1:length(protein.nms.picked.over2)){

			nms2 <- sapply(strsplit(protein.nms.picked.over2[kk], "___"), function(x){x[1]})

 
   			plot(xx2, yy2, main = paste(nms2, "peptides with", round(out.over2[kk,1],2), "\n (log2 scaled) with pval = ", signif(out.over2[kk, 6], 2)), ylab = nms.yy,   xlab = nms.xx)
   			abline(0,1, col = "grey")
			abline(cutoff, 1, col = "yellow")
			abline(-cutoff, 1, col = "yellow")



   			points(xx2[protein.id == protein.nms.picked.over2[kk]], yy2[protein.id == protein.nms.picked.over2[kk]], pch = 20, col = cols[kk])
		}
		dev.off()

	}

		


	# Let's make a table

	protein.id2 <- sapply(strsplit(rownames(res),  "___"), function(x){x[1]})
	protein.name <- sapply(strsplit(rownames(res),  "___"), function(x){x[2]})


	res2 <- data.frame("Accessions" = protein.id2, "Names" = protein.name, res)

	write.table (res2, paste(output.name, "_All.tsv", sep = ""), row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")


	



	if(!is.null(out) & nrow(out) > 0 ){

		protein.id2 <- sapply(strsplit(protein.nms.picked,  "___"), function(x){x[1]})
		protein.name <- sapply(strsplit(protein.nms.picked,  "___"), function(x){x[2]})


		out2 <- data.frame("Accessions" = protein.id2, "Names" = protein.name, out)

		write.table (out2, paste(output.name, "_Signif.tsv", sep = ""), row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")

		cat("done \n")


	}else{ cat("No protein is picked in this comparison \n") }


	

}



#########################################################################


do.iTRAQ()
