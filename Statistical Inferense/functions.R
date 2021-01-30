
pv <- function(X,theta){
  p <- numeric(nrow(X))
  for (i in 1:nrow(X)){ 
    x <- 1/(1+exp(-((X[i,]%*%theta))))
    p[i]<- x}
  return(p)}

pv_LOO <- function(X,theta){
  p <- numeric(999)
  for (i in 1:999){ 
    x <- 1/(1+exp(-((X[i,]%*%theta))))
    p[i]<- x}
  return(p)}

va <- function(p){
  v<- rep(0,length(p))
  for (i in 0:length(p)){ 
    v[i] <- p[i]*(1-p[i])}
  return(v)
}


D<- function(theta, y, X){diag(va(pv(X, theta)))}
L <- function(theta, dist,y, X){prod(dist(y, 1, 1/(1+exp(-(X%*%theta)))))}
l <- function(theta, y,size, X){log(L)}
S <- function(theta, y, X){t(X) %*% (y-pv(X,theta))}
I <- function(theta, y, X){t(X) %*% D(theta,y,X) %*% X}
NR <- function(theta, niter, y, X){
  for (i in 1:niter){
    theta<-theta+(solve((I(theta,y,X)))%*%S(theta,y,X))
    }
  return(theta)}

Std_er <- function(theta,y,X){sqrt(diag(solve(I(theta,y,X))))}
Se_theta.i <- function(theta,y,X,i){sqrt(diag(solve(I(theta,y,X)[i])))}


D_LOO<- function(theta, y, X){diag(va(pv_LOO(X, theta)))}
L_LOO <- function(theta, dist,y, X){prod(dist(y, 1, 1/(1+exp(-(X%*%theta)))))}
l_LOO <- function(theta, y,size, X){log(L)}
S_LOO <- function(theta, y, X){t(X) %*% (y-pv_LOO(X,theta))}
I_LOO <- function(theta, y, X){t(X) %*% D_LOO(theta,y,X) %*% X}


NR_LOO <- function(theta, niter, y, X){
  for (i in 1:niter){
    theta<-theta+(solve((I_LOO(theta,y,X)))%*%S_LOO(theta,y,X))
  }
  return(theta)}
