###############Distance and other small functions

#euclidean distance measure
eucldis <- function(x) sqrt(sum(x^2))

#L1 distance
L1_inf <- function(x) max(abs(x))

#Rescaling
rescale <- function(x){
  apply(x,MARGIN = 2,function(y){
    mn <- min(y);mx <- max(y);n <- length(y)
    a <- (n/(n-1)) * (mn - mx/n)
    b <- (n/(n-1)) * (mx - mn/n)
    (y - a)/(b - a)
  })
}

#### Non-Negative Function [,1:10]
noneg <- function(x){
  dtest <- !apply(exprs(x),1,function(row) any(row <= 0))
  exx <- exprs(x)[dtest==T,]
  paraa <- x@parameters
  describe <- x@description
  new("flowFrame",exprs=exx,parameters=paraa,description=describe)
}

##no na function[,1:11]
nona <- function(x){
  dtest <- !apply(exprs(x),1,function(row) any(is.na(row) | is.nan(row)))
  exx <- exprs(x)[dtest==T,]
  paraa <- x@parameters
  describe <- x@description
  new("flowFrame",exprs=exx,parameters=paraa,description=describe)
}

#log transformation
lnTrans <- function(x,notToTransform=c("SSC.W","TIME")){
  exx <- cbind(log(exprs(x)[,which(!(colnames(x) %in% notToTransform ))]), 
               exprs(x)[,which(colnames(x) %in% notToTransform)])
  colnames(exx) <- colnames(x)
  paraa <- x@parameters
  describe <- x@description
  new("flowFrame",exprs=exx,parameters=paraa,description=describe)
}