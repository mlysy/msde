library(msdeHeaders)
source("mOU-Kalman-revised.R")
source("r-code-parser.R")
require("mvtnorm")

# Only for mou case. This theta is one dimensional
# Need data.names to be X_1,...,X_ndims
test.post <- function(theta, Gamma, Lambda, Psi, init = 1,debug=FALSE){

	ndims <- length(Lambda) # length of Lambda == number of dimension

	# Initialize param.names
	param.names <- names(theta)

	# Initialize data.names
	data.names <- c("X1")
	for(nn in 2:ndims){
		data.names <- c(data.names,paste0('X',nn))
	}

	# Initialize the model file. Message will confirm the conversion
	quickParse(Gamma,Lambda,Psi)
	
	# Initialize model in C++
	message("Initializing Model in C++")
	moumod <- sde.make.model(ModelFile = "mOUModel.h",
													 param.names=param.names,
													 rebuild=TRUE)
	
	message("Model initialization complete.")
	nreps <- 0
	# check init
	if(is.numeric(init)){
		# Make a list of sde.init object
		# generate data
		nreps <- init
#		data.list <- NULL

		for(reps in 1:nreps){
			X0 <- runif(ndims, min = -1, max = 1)
			dT <- 1/250

			# generate nobs, burn
			burn <- 0
			nObs.sim <- 2000
			mou.sim <- sde.sim(model = moumod, x0 = X0, theta = theta,
												 dt = dT, dt.sim = dT, nobs = nObs.sim, burn=burn)

		# need to generate par.index,
			nObs.post <- 256
			nvar.obs <- c(ndims,sample(c(0:ndims),nObs.post-1,replace=TRUE))
			m <- sample(c(1:3),1)
			init.data <- sde.init(model = moumod,
														x = as.matrix(mou.sim$data[1:nObs.post,]),
														dt = mou.sim$dt,
														m = m,
														nvar.obs = nvar.obs,
														theta = theta)


		
			message("init.data completed.")

			nsamples <- 50000
			mou.post <- sde.post(model = moumod,
													 init = init.data,
													 hyper = NULL,
													 adapt = TRUE,
													 nsamples = nsamples,
													 burn = burn)

			h <- hist(mou.post$params, breaks=100, plot = FALSE)
			supp <- h$breaks
			res <- length(supp)

			ll.kalman <- rep(NA,res)
			len_WO <- dim(init.data$data)[1]

			data1 <- init.data$data[2:len_WO,]
			dts <- init.data$dt.m[2:(len_WO-1)]

			for(index in 1:res){
				TT <- supp[index]
				ll.kalman[index] <- dmou.Kalman(X = data1,
																				X0 = init.data$data[1,],
																				dt = dts,
																				dt0 = init.data$dt.m[1],
																				Gamma*TT,
																				Lambda*TT,
																				Psi*TT,
																				init.data$nvar.obs.m[2:len_WO],
																				log = TRUE,
																				debug=FALSE)
			}
			xden.kalman <- exp(ll.kalman - max(ll.kalman))
			h$density <- h$counts/sum(h$counts)
			h$density <- h$density/max(h$density)
			m.frac <- 1 - sum(init.data$nvar.obs.m)/(ndims*length(init.data$nvar.obs.m))
			main.title <- paste0("p",reps,": ndims=",ndims,", m.frac=",signif(m.frac,3))
			png(paste0('plot',reps,'.png'))
			plot(h,freq=FALSE,ylim=range(xden.kalman),xlim=range(supp),
						 xlab=expression(theta),ylab=expression(loglik(theta*" | X")),
						 main=main.title)
			lines(supp,xden.kalman,col="blue")
			abline(v=theta,lty=2,col='red')
			dev.off()
		}
	}
	else if(is.list(init)){
		stop("not implemented yet.")
	}
	else{
		stop("Please provide init as an integer or list of sde.init object.")
	}
	# define a variable to track the number of plots needed
#	return(mou.post)
}
