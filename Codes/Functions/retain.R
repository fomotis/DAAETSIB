# Decides if a file should be retiained or removed based on its status
# (from the good_measurement() function) as well as the number of particles
# it contains
#
# @param x dataframe result from good_measurement() function
# @param y dataframe result from good_measurement() function
# @param dil_x dilution level of x
# @param dil_y dilution level of y
#
# @return dataframe with status, dilution, filename of the retained file (either x or y)
#
# @examples
#
# \dontrun{
#   retain(goodFiles_200, goodFiles_500, "200", "500")
# }
#
# @export retain
#
retain <- function (x, y, dil_x = "100", dil_y = "200") {

  if(is.null(x) | is.null(y)) stop("Either of x or y is empty")
  #Dilution <- stringr::str_extract(as.character(bquote(x))[3], c(dil_x, dil_y))
  Dilution_x <- dil_x
  #Dilution <- stringr::str_extract(as.character(bquote(y))[3], c(dil_x, dil_y))
  Dilution_y <- dil_y
  needed <- data.frame()

  for( i in 1:nrow(x) ) {
    if(as.character(x$status[i]) == "good" & as.character(y$status[i]) == "good" ) {

      if(x$Size[i] > y$Size[i]) {

        needed <- rbind(needed, data.frame(x[i,], Dilution = Dilution_x ) )

      } else {

        needed <- rbind(needed, data.frame(y[i,], Dilution = Dilution_y ))
      }

    } else if(as.character(x$status[i]) == "good" & as.character(y$status[i]) == "bad") {

      needed <- rbind(needed, data.frame(x[i,], Dilution = Dilution_x) )

    } else if(as.character(x$status[i]) == "bad" & as.character(y$status[i]) == "good") {

      needed <- rbind(needed,data.frame(y[i,], Dilution = Dilution_y ))
    }
  }
  return(needed)
}
