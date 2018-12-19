# function to gate out bs4 from a debris filtered flowframe using a control flowframe.
# Called within the cell_debris function.
#
# @param flowframe supplied flowframe
# @param p1 same as cell_debris
# @param p2 same as cell_debris
# @param control_flowFrame control flowframe to used to identify debris
#
# @return list containing;
#     \itemize{
#         \item \strong{bs4bs5 -} debris filtered flowframe
#         \item \strong{others -} position of yet to be identified particles
#         \item \strong{debs -} position of debris
#         }
#
# @export debris_c


debris_c <- function(flowframe, p1, p2, control_flowFrame) {
  ##will gate out debris
  debris <- flowDensity::flowDensity(obj = flowframe, channels = c(p1, p2),
                                     position = c(F, NA), upper = c(T, T), use.upper = c(F, T),
                                     use.control = c(T, T), control = control_flowFrame)
  debs <- which(!is.na(debris@flow.frame@exprs[, 1])) #debris position
  others <- which(is.na(debris@flow.frame@exprs[, 1])) #other particles position

  #reduced_flowframe containing BS4 and BS5 only
  bs4bs5 <- flowframe[which(is.na(flowCore::exprs(debris@flow.frame)[, 1])), ]
  
  #plotting (plot1)
  flowDensity::plotDens(flowframe, c(p1, p2), xlab = p1, ylab = p2,
                        main = paste0("ID = ", flowCore::identifier(flowframe)), frame.plot = F)
  #Debris Filter
  points(debris@filter, type = "l", col = 2,lwd = 1,lty = 3)
  text(mean(debris@filter[, 1]), mean(debris@filter[, 2]), "Deb", col = 2)
  return(list(bs4bs5 = bs4bs5, debs = debs, others = others))

}
