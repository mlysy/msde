# Make the Heston SDE model
#mOU.model <- sde.make.model(ModelFile, PriorFile = NULL, data.names, params.names,...,debug=FALSE)

# The order in which the operations should be perform
#readLines(files)
#gsub()
#cat(c('d','o','g'),file="header.h", sep="")

# A is ndims by ndims matrix
# x is a ndims by 1 vector
# returns a string symbolic multiplication of Ax, ndims by 1 vector
matrixParse <- function(ndims, A, param = FALSE){
	astr <- matrix('', ncol = 1, nrow = ndims)
	# this loop is for row calculation
	for(nn in 1:ndims){
		# this loop is for column calculation of A
		for(coln in 1:ndims){
			astr[nn] <- paste0(astr[nn],signif(A[nn,coln],5),' * x[',(coln-1),']')
			if(param){
				astr[nn] <- paste0(astr[nn],' * theta[0]')
			}
			if((coln+1) <= ndims){
				astr[nn] <- paste0(astr[nn],' + ')
			}
		}
	}
	return(astr)
}

vectorParse <- function(ndims, avec, param=FALSE){
	astr <- matrix('',ncol = 1,nrow = ndims)

	for(nn in 1:ndims){
		astr[nn] <- paste0(astr[nn], signif(avec[nn],5))
		if(param){
			astr[nn] <- paste0(astr[nn],' * theta[0]')
		}
	}
	return(astr)
}

# This function parse the entire header file
# need to change: 
#		1. ndims
#   2. sdeDr - replace and/or increase/decrease the length of function
#		3. sdeDf - replace and/or increase/decrease the length of function
fileParse <- function(fileStr, strDr, strDf, ndims, debug=FALSE){
	# get the previous dimension
	dimStr <- grep('nDims',fileStr)
	oldDim <- as.numeric(gsub("[^0-9]","",fileStr[dimStr]))

	# replace previous dimension with current ndims
	fileStr[dimStr] <- gsub("[0-9]",ndims,fileStr[dimStr])
	if (debug) browser()
	# Process sdeDr
	tmpNew <- NULL
	rRowIndex <- grep("dr\\[",fileStr)
	tmpNew <- c(tmpNew,fileStr[1:(rRowIndex[1]-1)])

	oldhalf2 <- fileStr[(rRowIndex[length(rRowIndex)]+1):length(fileStr)]
	tmpNew <- c(tmpNew, strDr, oldhalf2)
	fileStr <- tmpNew
	# Done with sdeDr

	# Process sdeDf
	tmpNew <- NULL
	rRowIndex <- grep("df\\[",fileStr)
	tmpNew <- c(tmpNew,fileStr[1:(rRowIndex[1]-1)])

	oldhalf2 <- fileStr[(rRowIndex[length(rRowIndex)]+1):length(fileStr)]
	tmpNew <- c(tmpNew, strDf, oldhalf2)
	fileStr <- tmpNew
	# Done with sdeDf

	return(fileStr)
}


# this modifies the file
# void but modifies mOUModel.h
quickParse <- function(Gamma,Lambda,Psi){
	
	uPsi <- c(Psi) # vectorize the matrix
	ndims <- length(Lambda)
	
	# parse the calculation into string
	gamma.term <- matrixParse(ndims,Gamma,TRUE)
	lambda.term <- vectorParse(ndims,Lambda,TRUE)

	strDr <- matrix('',ncol=1,nrow=ndims)
	# sdeDr, the drift of SDE
	for(nn in 0:(ndims-1)){
		tmp <- paste0('  dr[',nn,'] = (',gamma.term[nn+1],' + ',lambda.term[nn+1],')')
		tmp <- paste0(tmp,';')
		strDr[nn+1] <- tmp
	}
	# sdeDf, the diffusion of SDE
	# Fill-in by column
	strDf <- matrix('',ncol=1,nrow=ndims*ndims)
	for(nn in 0:(ndims*ndims-1)){
		if(uPsi[nn+1] == 0){
			tmp <- paste0('  df[',nn,'] = ','0.0;')
		}
		else{
			tmp <- paste0('  df[',nn,'] = theta[0] * ',signif(uPsi[nn+1],5),';')
		}
		strDf[nn+1] <- tmp
	}

	headerFile <- readLines('headerT.h')
	tmp <- fileParse(headerFile, strDr, strDf, ndims)
	cat(tmp, file="mOUModel.h", sep='\n')
	message("Header file converted successfully.")
}
