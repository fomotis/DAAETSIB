
###TED Function
TED <- function(x, rescale = F, base_matrix = NULL, modification = NULL) {
  #rescale to 0-1 range
  if(rescale == T) x_rescale <- x else x_rescale <- rescale(x)
  orig_dist <- as.matrix(stats::dist(x_rescale, diag = F, upper = T))
  orig_dist <- as.numeric(orig_dist[base::upper.tri(orig_dist)])
  orig_density <- density(orig_dist)
  orig_trait <- orig_density$y
  
  if(is.null(base_matrix) & is.null(modification)) {
    
      stop("please supply a reference matrix or choose Modification accordingly")
    
  } else if(is.null(base_matrix) & modification = "M1") {
    
    dbar <- mean(orig_dist)
    vdbar <- (ncol(x_rescale) * (2 * nrow(x_rescale) + 3)) / 
      (90 * nrow(x_rescale) * (nrow(x_rescale) - 1))
    vd <- nrow(x_rescale) * vdbar
    b_trait <- dnorm(orig_density$x, mean = dbar, 
          sd = sqrt(vd))
    
  } else {
    
    base_dist <- as.matrix(stats::dist(base_matrix, diag = F, upper = T))
    base_dist <- as.numeric(base_dist[base::upper.tri(base_dist)])
    base_density <- density(base_dist)
    b_trait <- base_density$y
    
  }
  
  kldiv <- flexmix::KLdiv(cbind(orig_trait = orig_trait, b_trait = b_trait))[1, 2]
  ted <- 1 / (1 + kldiv)
  return(ted)
}