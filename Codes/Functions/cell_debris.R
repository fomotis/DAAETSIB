# Gates cyano cells of interest using 2 channels.
#
# @param flowframe flowframe to be filtered
# @param p1 character vector giving the first channel.
#           Defaults to RED.B.HLin (a measure of chlorophyll presence)
# @param p2 character vector giving the first channel.
#           Defaults to YEL.B.HLin (a measure of decaying pigments)
# @param positions dictates if both channels should be used in identifying the cell populations.
# @param control_use a logical flowframe indicating if a control flowframe be used for both channels.
#                    See \code{?flowDensity} for more
# @param control_flowFrame control flowframe
# @param interest character vector specifying if BS4, BS5 or both specie are of interest to be
#                 gated out of the sample
# @param mono_control is the control flowframe on first day if monocultures are being gated
# @param b4_control control flowframe for BS4 if a mix of BS4 and BS5 are being gated.
#                       Defaults to NULL
# @param b5_control control flowframe for BS5 if a mix of BS4 and BS5 are being gated.
#                       Defaults to NULL
# @return list containing the following;
# \itemize{
#         \item \strong{rframe -} reduced flowframe after debris have been filtered out.
#         \item \strong{fullframe -} full flowframe with indicator for debris in the expression matrix
#         \item \strong{Cell_count -} total number of cyano cells counted. This will be the number of BS4
#         if interest is BS4 and number of BS5 if interest is BS5 and addition of both
#         if interest is both.
#         \item \strong{N_particle -} is the total number of particles counted.
#         \item \strong{N_others -} is the total number of other particles counted.
#         }
#
#
#  @import Biobase
#  @export celldebris
#
#
celldebris <- function (flowframe, p1 = "RED.B.HLin", p2 = "YEL.B.HLin", positions = c(F, NA),
                        control_use = c(T, T), control_flowFrame = NULL,
                        interest = c("BS4", "BS5", "Both"),
                        mono_control = NULL, b4_control = NULL, b5_control = NULL) {

  if (interest != "Both" & is.null(mono_control)) {
    stop("Supply Control flowframe for mono_control")
  } else if((is.null(b4_control) | is.null(b5_control)) & interest == "Both"){
    stop("Supply Control flowframe for b4_control or b5_control")
  }

  #constructing the annotated data frame for the parameter
  if(interest == "Both"){
    ddata <- data.frame(rbind(flowframe@parameters@data,
                              c("BS4BS5.Indicator", "BS4BS5.Indicator", 1, 0, 1)))
    row.names(ddata) <- c(row.names(flowframe@parameters@data), '$P13')
  } else if(!is.null(mono_control)) {
    ddata <- data.frame(rbind(flowframe@parameters@data,
                              c("Debris.Indicator", "Debris.Indicator", 1, 0, 1)))
    row.names(ddata) <- c(row.names(flowframe@parameters@data), '$P13')

    #mono_control gating
    day1 <-  flowDensity::flowDensity(mono_control,
                                      channels = c(p1, p2),
                                      position = c(F, NA), upper = c(T, T),
                                      use.upper = c(F, T))
    day1_2 <- mono_control[which(is.na(flowCore::exprs(day1@flow.frame)[, 1])), ]
  } else stop("incorrect options supplied, check interest and mono_control")

  dvarMetadata <- flowframe@parameters@varMetadata
  ddimnames <- flowframe@parameters@dimLabels
  paraa <- Biobase::AnnotatedDataFrame(data = ddata, varMetadata = dvarMetadata,
                                       dimLabels = ddimnames)
  describe <- flowframe@description
  BS4BS5_ind <- rep(NA, flowCore::nrow(flowframe))

  if(!is.null(control_flowFrame)) {
    #gate debris
    debris <- debris_c(flowframe = flowframe, p1 = p1, p2 = p2,
                       control_flowFrame = control_flowFrame)
    bs4bs5 <- debris$bs4bs5
    others <- debris$others
    BS4BS5_ind[debris$debs] <- 0 #debris

    if(interest == "BS4"){
      #BS4b gating with control
      bs4_gate <- bs4_c(bs4bs5, p1, p2, day1_2, others)
      #reduced flow frame for BS4
      bs4_reduced <- bs4_gate$bs4_reduced
      #positions of the particles
      bs4_others <- bs4_gate$bs4_others
      others_nk <- bs4_gate$others_nk
      others_bs5 <- bs4_gate$others_bs5

      #RED.R.HLin and FSC.HLin gating
      #bs4_gate2 <- redfscgoodbad(bs4_reduced, p3, p4)
      #good <- bs4_gate2$good
      #bad <- bs4_gate2$bad
      #cut_off2 <- bs4_gate2$cut_off2

      #bs4_good <- bs4_others[good]
      #bs4_bad <- bs4_others[bad]
      #bs4_reduced2 <- bs4_reduced[good, ]
      #bs4_orig <- flowDensity::flowDensity(bs4_reduced2, c(p3, p4),
      #                      position = c(F, F), use.percentile = c(T, T),
      #                      percentile = c(0.95, 0.95), ellip.gate = T, scale = 0.99)
      #others_nk2 <- bs4_good[which(is.na(bs4_orig@flow.frame@exprs[, 1]))]
      #bs4_pos <- bs4_good[which(!is.na(bs4_orig@flow.frame@exprs[, 1]))]

      #plotting plot2
      #flowDensity::plotDens(bs4_reduced, c(p3, p4), xlab = p3,
      #                      ylab = p4, main = paste0("ID = ", flowCore::identifier(bs4_reduced)),
      #                      frame.plot = F)
      #abline(h = cut_off2, lty = 3, col = 3, lwd = 3)
      #points(bs4_orig@filter, type = "l", col = 2, lwd = 2, lty = 4)
      #text(mean(bs4_orig@filter[, 1]), mean(bs4_orig@filter[, 2]), "BS4", col = 2)

      #Indicators for each cell type
      BS4BS5_ind[others_bs5] <- 1 #others
      BS4BS5_ind[others_nk] <- 1 #others
      #BS4BS5_ind[others_nk2] <- 1 #others
      #BS4BS5_ind[bs4_bad] <- 1 #others
      BS4BS5_ind[bs4_others] <- 2 #BS4
      #number of BS4
      #n_BS4 <- flowCore::nrow(flowDensity::getflowFrame(bs4_reduced))
      #number of bs4s, bs5 and cyano
      n_debris <- sum(BS4BS5_ind == 0)
      n_others <- sum(BS4BS5_ind == 1)
      n_BS4 <- sum(BS4BS5_ind == 2)
      BS4_exx <- as.matrix(cbind(flowCore::exprs(flowframe),
                                 as.numeric(BS4BS5_ind)))
      #2 = BS4, 1 = others, 0 = debris
      colnames(BS4_exx)[length(colnames(BS4_exx))] <- 'Particle.Indicator'
      BS4_fflowframe <- new("flowFrame", exprs = BS4_exx, parameters = paraa, description = describe)
      #reduced flowframe
      reducedframe <- bs4_reduced#flowDensity::getflowFrame(bs4_orig)
      #full flowframe
      fflowframe <- BS4_fflowframe
      n_cyano <- n_BS4
    } else if(interest == "BS5"){
      bs5_gate <- bs5_c(bs4bs5, p1, p2, day1_2, others)
      #reduced flowframe for BS5
      bs5_reduced <- bs5_gate$bs5_reduced
      #positions of the particles
      others_nk <- bs5_gate$others_nk
      others_bs4 <- bs5_gate$others_bs4
      others_bs52 <- bs5_gate$others_bs52

      #RED.R.HLin and FSC.HLin gating
      #bs5_gate2 <- redfscgoodbad_bs5(bs5_reduced, p3, p4)
      #good <- bs5_gate2$good
      #bad <- bs5_gate2$bad
      #cut_off2 <- bs5_gate2$cut_off2

      #bs5_good <- others_bs52[good]
      #bs5_bad <- others_bs52[bad]
      #bs5_reduced2 <- bs5_reduced[good, ]
      #bs5_orig <- flowDensity::flowDensity(bs5_reduced2, c(p3, p4),
      #                                     position = c(F, F), use.percentile = c(T, T),
      #                                     percentile = c(0.95, 0.95), ellip.gate = T,
      #                                     scale = 0.99)
      #plotting plot2
      #flowDensity::plotDens(bs5_reduced, c(p3, p4), xlab = p3,
      #                      ylab = p4, main = paste0("ID = ", flowCore::identifier(bs5_reduced)),
      #                      frame.plot = F)
      #abline(h = cut_off2, lty = 3, lwd = 3, col = 2)
      #points(bs5_orig@filter, type = "l", col = 2, lwd = 2, lty = 4)
      #text(mean(bs5_orig@filter[,1]), mean(bs5_orig@filter[,2]), "BS5", col = 2)
      #n_BS5 <- flowCore::nrow(flowDensity::getflowFrame(bs5_reduced))

      #bs5 original positions
      #others_nk2 <- bs5_good[which(is.na(bs5_orig@flow.frame@exprs[, 1]))]
      #bs5_pos <- bs5_good[which(!is.na(bs5_orig@flow.frame@exprs[,1]))]

      #BS4, BS5 indicators
      BS4BS5_ind[others_bs4] <- 1 #others
      BS4BS5_ind[others_nk] <- 1 #others
      #BS4BS5_ind[others_nk2] <- 1 #others
      #BS4BS5_ind[bs5_bad] <- 1
      BS4BS5_ind[others_bs52] <- 2 #BS5
      #number of bs4s, bs5 and cyano
      n_debris <- sum(BS4BS5_ind == 0) #debris
      n_others <- sum(BS4BS5_ind == 1) #others
      n_BS5 <- sum(BS4BS5_ind == 2) #BS5
      #indicator for BS5
      BS5_exx <- as.matrix(cbind(flowCore::exprs(flowframe),
                                 as.numeric(BS4BS5_ind)))
      #2 = BS4, 1 = others,0 = debris
      colnames(BS5_exx)[length(colnames(BS5_exx))] <- 'Particle.Indicator'
      BS5_fflowframe <- new("flowFrame", exprs = BS5_exx, parameters = paraa, description = describe)
      #reduced flowframe
      reducedframe <- bs5_reduced#flowDensity::getflowFrame(bs5_orig)
      #full flowframe
      fflowframe <- BS5_fflowframe
      n_cyano <- n_BS5
    } else {
      #BS4
      bs4 <- flowDensity::flowDensity(bs4bs5, channels = c(p1, p2), position = c(NA, F),
                                      use.percentile = c(F, T),
                                      percentile = c(0.975, 0.975),use.control = c(T, T),
                                      control = c(b4_control, b4_control))
      #plotting
      points(bs4@filter, type = "l", lty = 4, col = 2, lwd = 1.5)
      text(mean(bs4@filter[,1]), mean(bs4@filter[,2]),"BS4",col = 2)

      bs4_reduced <- flowDensity::getflowFrame(bs4)
      n_BS4 <- flowCore::nrow(bs4_reduced)
      #bs5
      msg <- stringr::str_sub(capture.output(flowDensity::flowDensity(bs4bs5,
                                                                      channels = c(p1, p2),
                                                                      position=c(NA, F)))[1], 6, 27)
      if(stringr::str_detect(msg,"peak") == T |
         modes::bimodality_coefficient(flowCore::exprs(bs4bs5)[,p2]) < (5/9)){
        bs5 <- flowDensity::flowDensity(bs4bs5[which(is.na(flowCore::exprs(bs4@flow.frame)[, 1])),],
                                        channels = c(p1, p2),
                                        position = c(T, T), debris.gate = c(F, F),
                                        use.percentile = c(T, T), percentile = c(0.50,0.50),
                                        use.control = c(T, T),
                                        control = c(b5_control,b5_control))
      } else {
        bs5 <- flowDensity::flowDensity(bs4bs5[which(is.na(flowCore::exprs(bs4@flow.frame)[, 1])),],
                                        channels= c(p1, p2),
                                        position = c(T, T),
                                        use.percentile = c(F, T), percentile = c(0.50, 0.50),
                                        use.control=c(T, T),
                                        control=c(b5_control,b5_control))
      }

      #plotting
      points(bs5@filter,type = "l", lty = 4, col = 3)
      text(mean(bs5@filter[,1]), mean(bs5@filter[, 2]), "BS5", col = 2)

      bs5_reduced <- flowDensity::getflowFrame(bs5)
      n_BS5 <- flowCore::nrow(bs5_reduced)
      #indicator for bs4 and bs5
      bs4bs5_ind <- numeric(n_BS4+n_BS5)

      #reduced flowframe
      bs4bs5_ind[1:n_BS4] <- 1 #bs4=1, bs5=0
      #bind bs4 and bs5 together
      BS4BS5_rexx <- as.matrix(cbind(rbind(flowCore::exprs(bs4_reduced),
                                           flowCore::exprs(bs5_reduced)),bs4bs5_ind))
      colnames(BS4BS5_rexx)[length(colnames(BS4BS5_rexx))] <- 'BS4BS5.Indicator'
      BS4BS5_rframe <- new("flowFrame",exprs=BS4BS5_rexx,parameters=paraa,description=describe)
      reducedframe <- BS4BS5_rframe

      #full flowframe
      #BS4, BS5 indicators
      BS4BS5_ind <- NULL
      BS4BS5_ind[!is.na(flowCore::exprs(debris@flow.frame)[,1])] <- 0 #Debris
      ##BS4 and BS5 positions in the full flowframe
      bs4bs5_pos <- which(is.na(flowCore::exprs(debris@flow.frame)[,1]))
      BS4BS5_ind[bs4bs5_pos[!is.na(flowCore::exprs(bs5@flow.frame)[,1])]] <- 1 #BS5s
      bs4_pos <- which(is.na(flowCore::exprs(bs5@flow.frame)[,1]))
      BS4BS5_ind[bs4_pos[!is.na(flowCore::exprs(bs4@flow.frame)[,1])] ] <- 2 #BS4s

      BS4BS5_exx <- as.matrix(cbind(flowCore::exprs(flowframe),BS4BS5_ind ))
      colnames(BS4BS5_exx)[length(colnames(BS4BS5_exx))] <- 'Debris.Indicator'
      BS4BS5_fflowframe <- new("flowFrame",exprs=BS4BS5_exx,
                               parameters=paraa,description=describe)
      fflowframe <- BS4BS5_fflowframe
      n_cyano <- c(n_BS4, n_BS5)
    }

  } else stop("Supply control_flowframe")
  return(list(fullframe = fflowframe, rframe = reducedframe, Cell_count = n_cyano,
              N_particle = flowCore::nrow(flowframe), N_others = n_others))
}
