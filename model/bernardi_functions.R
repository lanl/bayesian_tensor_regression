##############
## Code developed by Bernardi 2023 to create the simulation data used
## in their experiments.
## Date added: April 2024
library(mnormt)

data.gen <- function(n,meanx,sdx,r){
  if(length(meanx)>1){
    grid = seq(0,1,length.out = r)
    X <- matrix(NA, r, n)
    for(i in 1:n) 
      X[,i] <- bs(grid, Boundary.knots = c(-0.1, 1.1), degree = 3,df = 15) %*% cbind(rmnorm(1, meanx, sdx))
  }
  if(length(meanx)==1) X = t(rnorm(n,meanx,sdx))
  return(t(X))
}

beta_conc <- function(xy) {
  out <- dbeta(.5+(xy[1]-xy[2])*3, 2,2)
}

beta_hist <- function(xy) ddirichlet(c(xy[2]-xy[1], 1-xy[2],xy[1]), alpha=c(4,4,4))

beta_tony <- function(xy,a,b,c,a2,b2,c2){
  x <- xy[1]
  y <- xy[2]
  res <- 0
  if(x<0.375 & y > .375){
    yy <- ((y-0.375)*1.6)
    xx <- x/.375
    res <-   ((xx)^(a-1)* yy^(b-1)*(1-(xx))^(b+c-1)*(1-yy)^(a+c-1))/((1-(xx)*yy)^(a+b+c))
  }
  if(x>0.375 & y < 0.75){
    xx <- ((x-0.375)*1.6)
    yy <- (y*4/3)
    res <- -(xx^(a2-1)*yy^(b2-1)*(1-xx)^(b2+c2-1)*(1-yy)^(a2+c2-1))/((1-xx*yy)^(a2+b2+c2))
  }
  res
}

beta_full <- function(xy) (0.3*(xy[2]-1/3))^3 + (0.3*(xy[1]-1/3))^3 - (0.2*(xy[2]-2/3))^2 - (0.2*(xy[1]-2/3))^2 

gen_beta <- function(method, s_seq, t_seq){
  
  st_seq  = expand.grid(s_seq, t_seq) 
  r1      = length(s_seq)
  r2      = length(t_seq)
  
  if(method == "quasi-concurrent"){
    out <- matrix(apply(st_seq,1, beta_conc), r1, r2)
    out <- out/sqrt(sum(out^2))
  }
  
  if(method == "historical"){
    out <- matrix(apply(st_seq,1, beta_hist), r1, r2)
    out <- out/sqrt(sum(out^2))
  }
  
  if(method == "local"){
    a1 <- 2
    b1 <- 13
    c1 <- 5
    a2 <- 13
    b2 <- 2
    c2 <- 5
    
    out <- matrix(apply(st_seq,1, beta_tony, a = a1, b = b1, c = c1, 
                        a2 = a2, b2 = b2, c2 = c2), r1, r2)
    out <- out/sqrt(sum(out^2))
  }
  if(method == "full"){
    out <- matrix(apply(st_seq,1, beta_full), r1, r2)
    out <- out/sqrt(sum(out^2))
  }
  
  return(out)
  
}

gen_alpha0 <- function(){
  alpha0     <- rnorm(15)
  alpha0     <- round(alpha0, 3)
  alpha0[4]  <- -0.2
  alpha0[15] <- 0
  alpha0[10] <- 0
  alpha0[11] <- 0.4
  alpha0[12] <- 0
  return(alpha0)
}

