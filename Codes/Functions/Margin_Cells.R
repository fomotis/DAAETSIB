# Removes or assign indicators to margin events
#
# @param flow.frame Flowframe containing margin events to be filtered out
# @param Channel The channel on which margin events are. Defaults to SSC.W (side scatter width)
# @param type The method to be used in gating out the margin cells. Can either be "manual" where user
#             supplies a cut off point on the channel, 1 = not margin 0 = margin
# @param  cut sould not be NULL if type = "manual"
#
# @return list containing;
# \itemize{
#         \item \strong{reducedflowframe -} flowframe without margin events
#         \item \strong{fullflowframe -} flowframe with an Margin.Indicator added as an
#         extra column added to the expression matrix to indicate which particles are margin events.
#         \item \strong{N_margin -} number of margin events recorded
#         \item \strong{N_cell -} numner of non-margin events
#         \item \strong{N_particle -} is the number of particles in total, i.e. N_cell + N_margin
#         }
# @examples
#
# \dontrun{
#    cellmargin(log_TransformedSet$Day1[[1]], "SSC.W", "estimate")
# }
# @importFrom methods new
# @export cellmargin

cellmargin <- function(flow.frame, Channel = "SSC.W", type = c("manual", "estimate"), cut = NULL){
  if(type=="manual" & !is.null(cut)) {
    margin.ind <- ifelse(flowCore::exprs(flow.frame)[,Channel] <= cut, T, F)
    n_margin <- sum(margin.ind == F)

    #plotting
    flowDensity::plotDens(flow.frame, c(Channel, "FSC.HLin"), main = flowCore::identifier(flow.frame))
    abline(v = cut, lwd = 1, lty = 4, col = 2)
    #constructing full flow frame with both margin and non-margin events
    exx <- as.matrix(cbind(flowCore::exprs(flow.frame),
                           as.numeric(margin.ind))) #1=not margin, 0=margin
    colnames(exx)[ncol(exx)] <- 'Margin.Indicator'
    #constructing the annotated data frame for the parameter
    ddata <- data.frame(rbind(flow.frame@parameters@data,
                              c("Margin.Indicator", "Margin Indicator", 1, 0, 1)))
    row.names(ddata) <- c(row.names(flow.frame@parameters@data), paste0("$", "P", ncol(exx)) )
    dvarMetadata <- flow.frame@parameters@varMetadata
    ddimnames <- flow.frame@parameters@dimLabels
    paraa <- Biobase::AnnotatedDataFrame(data = ddata, varMetadata = dvarMetadata,
                                         dimLabels = ddimnames)

    describe <- flow.frame@description
    fflowframe <- new("flowFrame", exprs = exx, parameters = paraa, description = describe)

    #constructing flowframe with only the non margin events
    exx2 <- exx[exx[,'Margin.Indicator'] == 1,]
    rflowframe <- new("flowFrame",exprs = exx2, parameters = paraa, description = describe)
  } else if(type == "estimate"){
    #bfilter <- density(exprs(flow.frame)[,Channel])
    #infl <- c(FALSE, diff(diff(bfilter$y) > 0) > 0)
    #if(length(unique(infl)) > 1){
    #  infl_point <- max(bfilter$x[infl],na.rm=T)
    #}else{
    #  infl_point <- max(bfilter$x,na.rm=T)
    #}
    infl_point1 <- flowDensity::deGate(flow.frame, Channel, bimodal = T)
    infl_point2 <- flowDensity::deGate(flow.frame, Channel, use.upper = T, upper = F)
    #plotting
    flowDensity::plotDens(flow.frame, c(Channel, "FSC.HLin"), main = flowCore::identifier(flow.frame))
    abline(v = infl_point1, lwd = 1, lty = 4, col = 2)
    abline(v = infl_point2, lwd = 1, lty = 4, col = 2)

    margin.ind <- ifelse(flowCore::exprs(flow.frame)[, Channel] > infl_point2 &
                           flowCore::exprs(flow.frame)[, Channel] < infl_point1,
                         T , F) #FALSE=Margin Event
    n_margin <- sum(margin.ind == F) #part of output

    #constructing full flow frame with both margin and non-margin events
    exx <- as.matrix(cbind(flowCore::exprs(flow.frame), as.numeric(margin.ind))) #1=not margin, 0=margin
    colnames(exx)[ncol(exx)] <- 'Margin.Indicator' #1 = not amrgin, 0 = margin

    #constructing the annotated data frame for the parameter
    ddata <- data.frame(rbind(flow.frame@parameters@data,
                              c("Margin.Indicator", "Margin Indicator", 1, 0, 1)))
    row.names(ddata) <- c(row.names(flow.frame@parameters@data), paste0("$", "P", ncol(exx)))
    dvarMetadata <- flow.frame@parameters@varMetadata
    ddimnames <- flow.frame@parameters@dimLabels
    paraa <- Biobase::AnnotatedDataFrame(data = ddata, varMetadata = dvarMetadata,
                                         dimLabels = ddimnames)
    describe <- flow.frame@description
    fflowframe <- new("flowFrame", exprs = exx, parameters = paraa, description = describe)

    #constructing flowframe with only the non margin events
    exx2 <- exx[exx[,'Margin.Indicator'] == 1, ]
    rflowframe <- new("flowFrame", exprs = exx2, parameters = paraa, description = describe)
  } else stop("Error: check your inputs")

  return(list(reducedflowframe = rflowframe, fullflowframe = fflowframe, N_margin = n_margin,
              N_cell = nrow(rflowframe), N_particle = nrow(fflowframe)))
}
