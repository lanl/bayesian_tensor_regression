##############
## Code to generate the B-spline toy simulation data.
## We mostly used this data to make sure all the dimensions were correct
## as we debugged the model.
## Authors: (mostly) JD Strait and (a bit from) AC Murph

library(plyr)
library(GIGrvg)
library(gtools)
library(coda)
library(fields)
library(splines)
library(glmnet)
library(LaplacesDemon)

################################################################################
## Generate simulated functional data (not distributional data to start)
## We can make this more complicated later, and use distributions, but I just did
## this to correctly match dimensions for the data structures we're working with in this model
################################################################################
# Recall that all the following must be defined here:
## n, K, txx, tyy, x, y, beta

# Initial settings
n = 200 # number of function pairs
K = 5 # number of discretized samples for each function

# Input/output grids (assume common dense grid to start)
txx = seq(from=0, to=1, length.out=K)
tyy = seq(from=0, to=1, length.out=K)

# Input functions
# Expanded from B-spline basis, equally-spaced knots, cubic polynomials
## M: Justin is using B-splines to simulate the 'input' functions.
##    Starting w getting the basis, then using random noise from a normal
##    for the weights.
degree   = 3
n_knots  = 15
bxx      = bs(txx, df=n_knots+degree, degree=degree)
nx_basis = dim(bxx)[2]
# matplot(txx, bxx, type="l") # basis functions

cxx   = matrix( rnorm(n*nx_basis, mean=0, sd=1), nrow=n, ncol=nx_basis )
# x_int = matrix( rnorm(n, mean=0, sd=0.5), nrow=n, ncol=K )
x     = cxx %*% t(bxx)
matplot(txx, t(x), type="l") # observed input functions

# True matrix-valued regression coefficient
# Use RBF kernel on txx, tyy
ls = 0.1
beta = matrix(NA, nrow=n, ncol=n)
for(i in 1:n){
  for(j in 1:n){
    # This assumes that effect of nearby quantiles far outweighs the effect
    # of quantiles further away.
    beta[i,j] = exp(-0.5*((txx[i]-tyy[j])^2) / (ls^2))
  }
}
image(beta) # image of plot of grid

# Output functions
sig = 0.3
eps = matrix( rnorm(K*n, mean=0, sd=sig), nrow=K, ncol=n )
beta0 = matrix(0, nrow=n, ncol=K)

# Closer quantiles are the ones with the effects, plus some noise.
y = beta0 + beta %*% x + t(eps)
# matplot(tyy, t(y), type="l") # observed output functions




