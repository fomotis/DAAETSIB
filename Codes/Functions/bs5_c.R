# function to gate out bs5 from a debris filtered flowframe using a control flowframe.
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
#         \item \strong{bs5_reduced -} reduced flowframe containing BS5
#         \item \strong{others_nk -} position of unknown particles
#         \item \strong{others_bs52 -} position of BS5s
#         \item \strong{others_bs4 -} position of BS4s, if present
# }
#
# @importFrom utils capture.output
# @export bs5_c
#
#

bs5_c <- function(bs4bs5, p1, p2, day1_2, others) {
  #BS4
  b4try <- try(flowDensity::flowDensity(bs4bs5, channels = c(p1, p2),
                                        use.upper = c(F, T), upper = c(T, F),
                                        position = c(NA, F)), silent = T)
  #reduced flow frame for BS4
  if(inherits(b4try, "try-error")){
    msg <- ""
    others_bs4 <- NULL
    others_bs5 <- others
    bs5_others <- bs4bs5
  } else {
    bs4 <- flowDensity::flowDensity(bs4bs5, channels = c(p1, p2), position = c(NA, F),
                                    use.upper = c(F, T), upper = c(T, F))
    others_bs4 <- others[which(!is.na(bs4@flow.frame@exprs[, 1]))]
    others_bs5 <- others[which(is.na(bs4@flow.frame@exprs[, 1]))]
    msg <- stringr::str_sub(capture.output(flowDensity::flowDensity(bs4bs5, channels = c(p1, p2),
                                                                    position = c(NA, F),
                                                                    use.upper = c(F, T),
                                                                    upper = c(T, F)))[1], 6, 27)
    bs5_others <- bs4bs5[which(is.na(bs4@flow.frame@exprs[, 1])), ]
  }

  #BS5
  if(stringr::str_detect(msg, "peak") == T |
     modes::bimodality_coefficient(flowCore::exprs(bs4bs5)[,p2]) < (5/9)){
    bs5s <- flowDensity::flowDensity(bs5_others, channels = c(p1, p2), position = c(T, T),
                                     use.percentile = c(T, T), percentile = c(0.02, 0.02),
                                     use.control = c(T, T), control = c(day1_2, day1_2),
                                     ellip.gate = T)
  } else {
    bs5s <- flowDensity::flowDensity(bs5_others, channels = c(p1, p2),
                                     position = c(T, T), use.percentile = c(F, T),
                                     percentile = c(0.10, 0.025), use.control = c(T, T),
                                     control = c(day1_2, day1_2),
                                     ellip.gate = T)
  }
  #potential bs5 alone
  bs5_reduced <- bs5_others[which(!is.na(bs5s@flow.frame@exprs[, 1])), ]
  others_nk <- others_bs5[which(is.na(bs5s@flow.frame@exprs[, 1]))]
  others_bs52 <- others_bs5[which(!is.na(bs5s@flow.frame@exprs[, 1]))]
  #plotting BS5
  points(bs5s@filter, type = "l", col = 2, lwd = 2, lty = 4)
  text(mean(bs5s@filter[,1]), mean(bs5s@filter[, 2]), "BS5", col = 2)

  return(list(bs5_reduced = bs5_reduced, others_nk = others_nk,
              others_bs4 = others_bs4, others_bs52 = others_bs52))
}
