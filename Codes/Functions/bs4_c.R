# function to gate out bs4 from a debris filtered flowframe using a control flowframe.
# Called within the cell_debris function.
#
# @param bs4bs5 debris filtered flowframe
# @param p1 same as cell_debris
# @param p2 same as cell_debris
# @param day1_2 debris gated mono_control
# @param others position of not yet identified particles
#
# @return
# \itemize{
#         \item \strong{bs4_reduced -} reduced flowframe containing BS4
#         \item \strong{others_nk -} position of unknown particles
#         \item \strong{bs4_others -} position of BS4s
#         \item \strong{others_bs5 -} position of BS5s, if present
#  }
#
#@export bs4_c

 bs4_c <- function(bs4bs5, p1, p2, day1_2, others) {
   bs5_try <- try(flowDensity::flowDensity(bs4bs5, channels = c(p1, p2),
                                           position = c(T, T),
                                           use.percentile = c(T, T),
                                           percentile = c(0.95, 0.95)), silent = T)
   if(inherits(bs5_try, "try-error")){
     others_bs5 <- NULL
     others_bs4 <- others
     bs4s <- bs4bs5
   } else {
     bs5 <- flowDensity::flowDensity(bs4bs5, channels = c(p1, p2),
                                     position = c(T, T),
                                     use.percentile = c(T, T),
                                     percentile = c(0.95, 0.95))
     others_bs5 <- others[which(!is.na(bs5@flow.frame@exprs[, 1]))]
     others_bs4 <- others[which(is.na(bs5@flow.frame@exprs[, 1]))]
     bs4s <- bs4bs5[which(is.na(flowCore::exprs(bs5@flow.frame)[, 1])), ]
   }
   #BS4
   bs4_1 <- flowDensity::flowDensity(bs4s, channels = c(p1, p2), position = c(NA, F),
                                     use.percentile = c(F, T),
                                     percentile = c(0.975, 0.975), use.control = c(T, T),
                                     control = c(day1_2, day1_2), ellip.gate = T)
   others_nk <- others_bs4[which(is.na(bs4_1@flow.frame@exprs[, 1]))]
   bs4_others <- others_bs4[which(!is.na(bs4_1@flow.frame@exprs[, 1]))]
   #plotting
   points(bs4_1@filter, type = "l", col = 2, lwd = 2, lty = 4)
   text(mean(bs4_1@filter[, 1]), mean(bs4_1@filter[, 2]), "BS4", col = 2)

   #reduced flowframe for BS4
   bs4_reduced <- bs4s[which(!is.na(bs4_1@flow.frame@exprs[, 1])), ]
   return(list(bs4_reduced = bs4_reduced, others_nk = others_nk, bs4_others = bs4_others,
               others_bs5 = others_bs5))
 }
