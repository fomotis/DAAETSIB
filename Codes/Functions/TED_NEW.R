
#source("distances.R")

base_traitmatrix_new <- function(x, rescale = F, n_gen = 1500){
  
  #rescale to 0-1 range
  if(rescale == F) x_rescale <- x else x_rescale <- rescale(x)
  
  #compute convex hull and depth
  cx_A <- geometry::convhulln(x_rescale, "FA")
  #x_meds <- mrfDepth::hdepthmedian(x_rescale, maxdir = 3000)$median
  #points of convex hull
  #cx_points <- cx_A$hull
  
  #filter out observation that makes the convex hull
  #hull_points <- unique(x_rescale[t(cx_points), ])
  
  #compute minimum and maximum distance of points on the convex hull to the depth-median vector
  #hull_distance <- apply(hull_points - x_meds, 1, eucldis)
  #minmax <- c(min(hull_distance), max(hull_distance))
  
  N <- 0
  test_data <- data.frame()
  #mn <- apply(x_rescale, 2, min); mx <- apply(x_rescale, 2, max)
  while(N <= nrow(x_rescale)){
    #generate data from uniform distribution
    test <- runif(ncol(x_rescale), apply(x_rescale, 2, min),
          apply(x_rescale, 2, max))
    
    #compute their distance from the inner depth
    #test_distance <- eucldis(test - x_meds)
    
    #compute the area of the convex-hull if they are added to the data
    test_A <- geometry::convhulln(rbind(x_rescale, test), "FA")$area
    if(test_A == cx_A$area) {
      test_data <- rbind(test_data, test)
      N <- N + 1
    } else next
    if(N == n_gen) break
  }
  
  return(testdataS)
}