########################################################################################
# Original data must be a .csv file with columns: id, day, antigen, rep1, rep2, ...    #
########################################################################################

########################################################################################
#### Functions to implement DFR(eq) method							   #
#### Last update: March 17, 2025                                                        #
########################################################################################

library(dplyr)

######################
# Helper functions:  #
######################

impose.mono <- function(x, dir){
	# Function to impose monotonicity on elements of the vector, x
	# dir="incr" gives increasing monotonicity
	# dir="decr" gives decreasing monotonicity
  
	if(dir=="decr"){ x <- rev(x) }
	for(i in 2:length(x)){
		idx <- is.na(x[(i-1):i])
		if (sum(idx)==1){
			x[i] <- x[(i-1):i][!idx]
		} else {
			x[i] <- max(x[i-1],x[i])
		}
	}
	if(dir=="decr"){ x <- rev(x) }
      x
}

perm <- function(dat, nExp, nCtl){
	N <- nExp+nCtl 	  # total number of exp and neg ctl wells
	k <- NROW(dat)      # number of peptide pools
	B <- choose(N,nExp) # number of perms needed for complete enumeration

	if(B < 20) stop("Too few replicates to use this method (B < 20).")

	mu.e <- rowMeans(dat[,1:nExp], na.rm=TRUE)        # vector of peptide pool means
	mu.c <- rowMeans(dat[,nExp+(1:nCtl)], na.rm=TRUE) # vector of neg ctl means

	# test statistic for observed data
	t <- mu.e-mu.c                  
	t.sort <- sort(t)
	index <- order(t)

	# samp matrix contains all possible permutations of columns of dat matrix:
	samp <- expand.grid(rep(list(0:1),N))
	samp <- samp[apply(samp,1,sum)==nExp,]
	samp <- as.matrix(samp)

	tPerm <- matrix(0,nrow=k,ncol=B)
  
	# calculate test statistics for each perm sample in samp:
	for(i in 1:B){
		perm.dat <- dat[,order(samp[i,])]
		mu.exp <- rowMeans(perm.dat[,(1:nExp)+nCtl], na.rm=TRUE)
		mu.ctl <- rowMeans(perm.dat[,1:nCtl], na.rm=TRUE)
		tPerm[,i] <- mu.exp-mu.ctl
	}
  
	# order rows of tPerm to correspond with sorted test statistics in t.sort:
	if(k==1){ tPerm.mono <- tPerm.sort <- tPerm }
	if(k>1){
		tPerm.sort <- tPerm[index,]    
		tPerm.mono <- apply(tPerm.sort,2,impose.mono,dir="incr")
	}
  
	# calculate adjusted p-values:
	tpvalue <- apply(tPerm.mono>t.sort | near(tPerm.mono,t.sort) , 1, mean, na.rm=TRUE)
	
	# enforce monotonicity on adjusted p-values in tpvalue:
	if(k>1){ tpvalue <- impose.mono(tpvalue, dir="decr") }
        
	list(tstat=t, tadjp=tpvalue[order(index)])
}

#####################
# Main function:    #
#####################

elsdfreq <- function(data, nExp, nCtl, nameCtl, alpha=0.05){
	data <- as.data.frame(data)
	data <- data[order(data[,1],data[,2]),]
	ncdat <- data[data[,3]==nameCtl,]
	dat <- data.frame(matrix(NA,nrow=NROW(data[data[,3]!=nameCtl,]),ncol=3+nExp+nCtl))
	dat[,1:(3+nExp)] <- data[data[,3]!=nameCtl,]
	names(dat) <- c(names(data)[1:(3+nExp)],paste("c",1:nCtl,sep=""))
	for(id in unique(data[,1])){
	for(day in unique(data[,2])){
		dat[dat[,1]==id & dat[,2]==day,3+nExp+(1:nCtl)] <- ncdat[ncdat[,1]==id & ncdat[,2]==day,-(1:3)]
	}}
	
	### dat = a data.frame with ptid, day, antigen in cols 1 to 3,
	### peptide replicates in col4:(3+nExp), control replicates in
	### col(4+nExp):col(3+nExp+nCtl)
	### Same control replicates are repeated for each peptide per unique ptid

	### nExp = no. of experimental wells per peptide
	### nCtl = no. of negative control wells
	### Use alpha=0.05 as default
  
	ind <- unique(dat[,1:2])
	adjp <- teststat <- pos.call <- NULL

	# perform permutation test on each unique ptid*day combination with
	# multiplicity adjustment for the number of peptide pools
  
	for(i in 1:NROW(ind)){
		temp.dat <- as.matrix(dat[dat[,1]==ind[i,1] & dat[,2]==ind[i,2],-(1:3)])
		temp.id <- ind[i,1]
		temp.day <- ind[i,2]
		temp.pep <- dat[dat[,1]==ind[i,1] & dat[,2]==ind[i,2],3]
		if (length(temp.pep)==1) temp.dat <- t(temp.dat)
		temp.adjp <- temp.teststat <- rep(NA, NROW(temp.dat))
 
		nas <- is.na(temp.dat)
		nCtlc <- nCtl - sum(nas[1,nExp+(1:nCtl)]) # observed number of neg ctl reps
		nExpC <- nExp - apply(nas[,1:nExp],1,sum) # observed numbers of exp reps
		idx <- (nExpC>=3 & nCtlc>=3) | (nExpC>=2 & nCtlc>=4)
		if (sum(idx)>0){
			temp.dat <- as.matrix(temp.dat[idx,])
			if (sum(idx)==1 & length(idx)>1) temp.dat <- t(temp.dat)
			out <- perm(dat=temp.dat, nExp=nExp, nCtl=nCtl)
			temp.teststat[idx] <- out$tstat
			temp.adjp[idx] <- out$tadjp
		}
		teststat <- c(teststat, temp.teststat)
		adjp <- c(adjp, temp.adjp)
				
      	if (nCtlc<nCtl) warning(paste("ptid", temp.id, "has", nCtl-nCtlc, "missing value(s) for negative controls", "on day", temp.day))
		if (min(nExpC)<nExp){
			for (j in 1:sum(nExpC<nExp)){
				warning(paste("ptid", temp.id, "has", nExp-nExpC[nExpC<nExp][j], "missing value(s) for", temp.pep[nExpC<nExp][j], "on day", temp.day))
			}
		}
	}
	pos.call <- ifelse(adjp <= alpha,1,0)    

	output <- data.frame(cbind(dat,round(teststat,5),round(adjp,5),pos.call))
	colnames(output)[(ncol(dat)+1):(ncol(dat)+3)] <- c("t-stat","adjp","pos")
	output
}

