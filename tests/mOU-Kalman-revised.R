# The full model dYt = (Gamma * Yt + Lambda) dt + Phi dBt
# Yt = (Xt, Zt), where Xt's are observed and Zt's are latent
# With the above parameterization and discretization, we get that
# Y_{n+1} | Yn ~ N(Yn + (Gamma * Yn + Lambda) dt, Phi Phi' dt)

# functions that parameterize Gamma, Lambda, Phi
# Gamma, Lambda, Phi are known functions of unknown theta

# Kalman filter method

# Calculates the conditional mean and variance for Zn | X_{0:N}.
# This function is only called when par.index is not ndim, some parts of Y are latent.
# @param X is the n-th data point
# @param par.index is a scalar value
# @param mu and sigma are the mean and variance of Yn|Y0
mOU.MV <- function(X, mu, sigma, par.index, ndim, debug = FALSE){
	qn <- par.index # qn
	# debugger check
	if(debug) browser()
	X_obs <- X[1:qn] # selecting the observed Xn's
	m2 <- mu[1:qn]
	m1 <- mu[(qn+1):ndim]
	# double check the elements selection
	s22 <- sigma[1:qn,1:qn] # qn x qn
	s21 <- sigma[1:qn,(qn+1):ndim] # qn x (d - qn)
	s12 <- sigma[(qn+1):ndim,1:qn] # (d - qn) x qn
	s11 <- sigma[(qn+1):ndim,(qn+1):ndim] # (d - qn) x (d - qn)
	# possible for efficiency improvement
	s22inv <- solve(s22,diag(qn)) # solving the rref with the right dimension (i think)
	m_n <- m1 + s12 %*% s22inv %*% (X_obs - m2) # formula calculation for m_n and V_n
	V_n <- s11 - s12 %*% s22inv %*% s21
	MV <- list(M = m_n, V = V_n, obsM = m2, obsV = s22)
	return(MV)
}


# PDF computation through Kalman filter
# @param X \code{N x d} matrix of partially observed data at intervals of \code{dt}. This implies the i-th row is the i-th observation.
# @param dt is either a scalar or a vector for multiresolution
# @param Gamma is a \code{d x d} matrix
# @param Lambda is a \code{d}-dimensional vector
# @param Phi is a \code{d x d} matrix
# @param par.index indicates which of the corresponding row in X are observed and how many are observed.

dmou.Kalman <- function(X, X0, dt, dt0, Gamma, Lambda, Phi, par.index, log = TRUE, debug = FALSE){
	# dimension calculations
	ndim <- ncol(X)
	nObs <- nrow(X)

	# pre-allocate space for p(X)
	X.pdf <- 0
	# pre-allocate memory for dts (multiresolution purpose)
	dts <- rep(NA,nObs)

	# define the resolution for calculations
	if(length(dt) == 1){
		dts <- rep(dt, nObs)
	}
	else	if(length(dt) == (nObs-1)){
		dts <- dt
	}
	else{
		stop("Incorrect dt.")
	}

	# initialize the A_n, b_n, C_n
	A_n <- diag(ndim) + Gamma * dt0
	b_n <- Lambda * dt0
	C_n <- Phi * sqrt(dt0)

	# mu and Sigma calculation
	# Y1|Y0 = (X1,Z1) | Y0 ~ N(mu_1, Sigma_1)
	mu_n <- A_n %*% X0 + b_n
	Sigma_n <- crossprod(C_n)
	
	for(nn in 1:nObs){
	# coming into the loop we already calculated mu_n and Sigma_n which are needed for p(Xn | X_{0:n-1})
		qn <- par.index[nn]
	# conditional m and V calculation for Z_n, ie Zn | Xn,X0 ~ N(m_n,V_n)
		if (qn == ndim){

			X_obs <- X[nn,]
			
			if(qn == 1) Sigma_n <- matrix(Sigma_n,ncol=1,nrow=1)
#			X.pdf[nn-1] <- dmvnorm(X_obs, mean = mu_n, sigma = t(Sigma_n), log = log)
			X.pdf <- X.pdf + dmvnorm(X_obs, mean = mu_n, sigma = t(Sigma_n), log = log)

			# mu and Sigma calculation
			mu_n <- A_n %*% X_obs + b_n # mu_n is now mu_{n+1}
			Sigma_n <- crossprod(C_n) # V_n is now V_{n+1}
		}
		else if (qn != 0){

			# calculates the conditional mean and variance, m_n and V_n throught the function
			X_obs <- X[nn,1:qn]

			condMV <- mOU.MV(X[nn,], mu_n, Sigma_n, qn, ndim, debug = FALSE)
			m_n <- condMV$M
			V_n <- condMV$V
			m1 <- condMV$obsM
			s11 <- condMV$obsV

		# the density for p(Xn | X_{0:n-1}) using qn entries of mu_n and Sigma_n
			X.pdf <- X.pdf + dmvnorm(X_obs, mean = mu_n[1:qn], sigma = matrix(Sigma_n[1:qn,1:qn],qn,qn), log=log)
		# update mu_n and Sigma_n using m_n and V_n
			mu_n <- b_n + A_n %*% c(X_obs,m_n)# mu_n is now mu_{n+1}
			Sigma_n <- crossprod(C_n) + A_n[,(qn+1):ndim] %*% V_n %*% t(A_n[,(qn+1):ndim]) # V_n is now V_{n+1}
		}
		else{ # for the case when qn = 0 
			m_n <- mu_n
			V_n <- Sigma_n
			# progress the mu and sigma
			mu_n <- A_n %*% m_n + b_n
			Sigma_n <- crossprod(C_n) + A_n %*% V_n %*% t(A_n)
		}

		A_n <- diag(ndim) + Gamma * dts[nn]
		b_n <- Lambda * dts[nn]
		C_n <- Phi*sqrt(dts[nn])
	}
	if(debug) browser()
	# check for log = TRUE
	ll <- X.pdf
	if(!log) ll <- exp(ll)
	return(ll)
}

