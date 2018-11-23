
###TED Function
TED <- function(x, rescale = F, base_matrix = NULL){
  if(is.null(base_matrix)) stop("please supply a reference matrix")
  
  base_dist <- as.matrix(stats::dist(base_matrix, diag = F, upper = T))
  base_dist <- as.numeric(base_dist[base::upper.tri(base_dist)])
  
  #rescale to 0-1 range
  if(rescale == T) x_rescale <- x else x_rescale <- rescale(x)
  
  orig_dist <- as.matrix(stats::dist(x_rescale, diag = F, upper = T))
  orig_dist <- as.numeric(orig_dist[base::upper.tri(orig_dist)])
  
  kldiv <- flexmix::KLdiv(cbind(b_trait = density(base_dist)$y,
                                densityd_trait = density(orig_dist)$y ))[1, 2]
  ted <- 1 / (1 + kldiv)
  return(ted)
}