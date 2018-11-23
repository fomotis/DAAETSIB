TOP_new <- function(trait_matrix, n_gen = nrow(trait_matrix)){
  #compute convex hull and depth
  cx_A <- geometry::convhulln(trait_matrix, "FA")
  x_meds <- mrfDepth::hdepthmedian(trait_matrix, maxdir = 3000)$median
  
  #points of convex hull
  cx_points <- cx_A$hull
  
  #filter out observation that makes the convex hull
  #hull_points <- unique(trait_matrix[t(cx_points), ])
  
  #compute minimum and maximum distance of points on the convex hull to the depth-median vector
  #hull_distance <- apply(hull_points-x_meds, 1, eucldis)
  #minmax <- c(min(hull_distance), max(hull_distance))
  
  N<-0
  test_matrix <- data.frame()
  while(N <= nrow(trait_matrix)){
    #generate data from uniform distribution
    test <- runif(ncol(trait_matrix), apply(trait_matrix, 2, min),
                  apply(trait_matrix, 2, max))
    
    #compute their distance from the inner depth
    #test_distance <- eucldis(test - x_meds)
    
    #compute the area of the convex-hull if they are added to the data
    test_A <- geometry::convhulln(rbind(trait_matrix,test), "FA")$area
    if(test_A == cx_A$area){
      test_matrix <- rbind(test_matrix, test)
      N <- N + 1
    } else next
    if(N == n_gen) break
  }
  #calculate the total number of convex hull that can be peeled off
  N_possib <- 0
  while(is.null(test_matrix) == F & (nrow(test_matrix) > dim(test_matrix)[2])){
    mat_reduce <- t(geometry::convhulln(test_matrix, "Pp"))
    test_matrix <- unique(test_matrix[-mat_reduce,])
    if(is.null(nrow(test_matrix)) == T) break
    N_possib <- N_possib + 1
  }
  
  #calculating the classical TOP
  area <- c(); row_num <- c()
  while(is.null(trait_matrix)==F & (nrow(trait_matrix) > dim(trait_matrix)[2])){
    hulls <- geometry::convhulln(trait_matrix,"FA")
    area <- c(area,hulls$area)
    mat_reduce <- t(hulls$hull)
    row_num <- c(row_num, nrow(unique(trait_matrix[mat_reduce,])))
    trait_matrix <- unique(trait_matrix[-mat_reduce, ])
    if(is.null(nrow(trait_matrix)) == T) break
  }
  return(list(TA = sum(area)/N_possib, Areas = area, N = N_possib,  N_points = row_num, 
              TA2 = sum((area/N_possib) * row_num)))
}