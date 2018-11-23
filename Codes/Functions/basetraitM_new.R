
#source("distances.R")

base_traitmatrix_new <- function(x, rescale = F, n_gen = 1500){
  
  #rescale to 0-1 range
  if(rescale == F) x_rescale <- x else x_rescale <- rescale(x)
  
  #compute area of the convex hull
  cx_A <- geometry::convhulln(x_rescale, "FA")
  N <- 0
  test_data <- data.frame()
  while(N <= nrow(x_rescale)){
    #generate data from uniform distribution
    test <- runif(ncol(x_rescale), apply(x_rescale, 2, min),
          apply(x_rescale, 2, max))
    
    #compute the area of the convex-hull if they are added to the data
    test_A <- geometry::convhulln(rbind(x_rescale, test), "FA")$area
    if(test_A == cx_A$area) {
      test_data <- rbind(test_data, test)
      N <- N + 1
    } else next
    if(N == n_gen) break
  }
  
  return(test_data)
}