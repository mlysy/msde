## ---- eval = FALSE-------------------------------------------------------
#  sde.drift <- function(x, theta) {
#    dr <- c(theta[1]*x[1] - theta[2]*x[1]*x[2], # alpha * H - beta * H*L
#        theta[2]*x[1]*x[2] - theta[3]*x[2]) # beta * H*L - gamma * L
#    dr
#  }

## ---- eval = FALSE-------------------------------------------------------
#  sde.diff <- function(x, theta) {
#    df <- matrix(NA, 2, 2)
#    df[1,1] <- theta[1]*x[1] + theta[2]*x[1]*x[2] # alpha * H + beta * H*L
#    df[1,2] <- -theta[2]*x[1]*x[2] # -beta * H*L
#    df[2,1] <- df[1,2] # -beta * H*L
#    df[2,2] <- theta[2]*x[1]*x[2] + theta[3]*x[2] # beta * H*L + gamma * L
#    df
#  }

