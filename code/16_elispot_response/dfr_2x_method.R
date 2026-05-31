########################################################################################
# Original data must be a .csv file with columns: id, day, antigen, rep1, rep2, ...    #
########################################################################################

########################################################################################
#### Functions to implement DFR(2x) method							   #
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

# functions to calculate means of resampled exp and nCtl data
bf <- function(x,dat,n,B){
	row <- na.omit(as.numeric(dat[x,1:n]))
	rowMeans(matrix(sample(row,size=length(row)*B,replace=TRUE),nrow=B))
}

bcf <- function(x,dat,ne,nc,B){
	row <- na.omit(as.numeric(dat[x,ne+(1:nc)]))
	rowMeans(matrix(sample(row,size=length(row)*B,replace=TRUE),nrow=B))
}

bootfn <- function(dat, nExp, nCtl, ho.c, B, seed){
	k <- NROW(dat)    # number of peptide pools

	mu.e <- rowMeans(dat[,1:nExp], na.rm=TRUE)        # vector of peptide pool means
	mu.c <- rowMeans(dat[,nExp+(1:nCtl)], na.rm=TRUE) # vector of neg ctl means

	# test statistic for observed data
	t <- mu.e - mu.c - ho.c
	t.sort <- sort(t)
	index <- order(t)

	set.seed(seed)

	# Calculate bootstrap means for exp and nCtl
	if (k==1){
		rowe <- na.omit(as.numeric(dat[,1:nExp]))
		rowc <- na.omit(as.numeric(dat[,nExp+(1:nCtl)]))
		bemeans <- rowMeans(matrix(sample(rowe,size=length(rowe)*B,replace=TRUE),nrow=B))
		bcmeans <- rowMeans(matrix(sample(rowc,size=length(rowc)*B,replace=TRUE),nrow=B))
	} else {
		bemeans <- t(apply(matrix(1:k),1,bf,dat,nExp,B))
		bcmeans <- t(apply(matrix(1:k),1,bcf,dat,nExp,nCtl,B))
	}
	tBoot <- (bemeans - mu.e) - (bcmeans - mu.c)

	# order rows of tBoot to correspond with sorted test statistics in t.sort:
	if (k>1){
		tBoot.sort <- tBoot[index,]
		# enforce monotonicity on columns of tBoot.sort:
            tBoot.mono <- apply(tBoot.sort,2,impose.mono,dir="incr")
	}

	# calculate adjusted p-values:
	if(k==1){
		tpvalue <- mean(tBoot > t.sort | near(tBoot,t.sort))
	} else {
	  tpvalue <- rowMeans(tBoot.mono > t.sort | near(tBoot.mono,t.sort))
	}

	# enforce monotonicity on adjusted p-values in tpvalue:
	if (k>1){ tpvalue <- impose.mono(tpvalue,dir="decr") }

	list(tstat=t, tadjp=tpvalue[order(index)])
}

#####################
# Main function:    #
#####################

elsdfr2x <- function(data, nExp, nCtl, nameCtl, ho.c=log10(2), alpha=0.05, B=1000, seed=9456845){
	data <- as.data.frame(data)
	data <- data[order(data[,1],data[,2]),]
	ncdat <- data[data[,3]==nameCtl,]
	dat <- data.frame(matrix(NA,nrow=NROW(data[data[,3]!=nameCtl,]),ncol=3+nExp+nCtl))
	dat[,1:(3+nExp)] <- data[data[,3]!=nameCtl,1:(3+nExp)]
	names(dat) <- c(names(data)[1:(3+nExp)],paste("c",1:nCtl,sep=""))
	for(id in unique(data[,1])){
	for(day in unique(data[,2])){
		dat[dat[,1]==id & dat[,2]==day,3+nExp+(1:nCtl)] <- ncdat[ncdat[,1]==id & ncdat[,2]==day,-(1:3)]
	}}
	idx <- which(rowSums(dat[,-(1:3)]==0, na.rm=TRUE) > 0)
	if (length(idx)>0){ dat[idx,-(1:3)] <- dat[idx,-(1:3)] + 1 }
	dat[,-(1:3)] <- log10(dat[,-(1:3)])

	### dat = a data.frame with ptid, day, antigen in cols 1 to 3,
	### peptide replicates in col4:(3+nExp), control replicates in
	### col(4+nExp):col(3+nExp+nCtl)
	### Same control replicates are repeated for each peptide per unique ptid*day

	### nExp = no. of experimental wells per peptide
	### nCtl = no. of negative control wells
	### ho.c = value for null hyp: Ho: mu.e - mu.c = ho.c
	### Use alpha=0.05 as default
	### B = no. of bootstrap samples, use 1000 as default

	ind <- unique(dat[,1:2])
	adjp <- teststat <- pos.call <- NULL

	# perform bootstrap test on each unique ptid*day combination with
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
		nExpC <- nExp - apply(nas[,1:nExp],1,sum,na.rm=TRUE) # observed numbers of exp reps
		idx <- (nExpC>=3 & nCtlc>=3) | (nExpC>=2 & nCtlc>=4)
		if (sum(idx)>0){
			temp.dat <- as.matrix(temp.dat[idx,])
			if (sum(idx)==1 & length(idx)>1) temp.dat <- t(temp.dat)
			out <- bootfn(dat=temp.dat, nExp=nExp ,nCtl=nCtl, ho.c=ho.c, B=B, seed=seed)
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

	output <- data.frame(cbind(dat[,1:3],round(teststat,5),round(adjp,5),pos.call))
	colnames(output) <- c("ptid","day","peptide","t-stat","adjp","pos")
	output
}
