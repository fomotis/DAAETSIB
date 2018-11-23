########function to screen out margin events
#### Flagging Margin Events
##1 = not margin 0=margin

cellmargin <- function(flow.frame,Channel="SSC.W",type=c("manual","estimate"),cut=NULL){
  if(type=="manual" & !is.null(cut)){
    margin.ind <- ifelse(exprs(flow.frame)[,Channel] <= cut, T, F)
    n_margin <- sum(margin.ind == F)
    
    #constructing full flow frame with both margin and non-margin events
    exx <- as.matrix(cbind(exprs(flow.frame),as.numeric(margin.ind))) #1=not margin, 0=margin
    colnames(exx)[ncol(exx)] <- 'Margin.Indicator'
    #constructing the annotated data frame for the parameter
    ddata <- data.frame(rbind(flow.frame@parameters@data,
                              c("Margin.Indicator","Margin Indicator",1,0,1)))
    row.names(ddata) <- c(row.names(flow.frame@parameters@data),paste0("$","P",ncol(exx)) )
    dvarMetadata <- flow.frame@parameters@varMetadata
    ddimnames <- flow.frame@parameters@dimLabels
    paraa <- Biobase::AnnotatedDataFrame(data=ddata,varMetadata=dvarMetadata,dimLabels=ddimnames)
    
    describe <- flow.frame@description
    fflowframe <- new("flowFrame",exprs=exx,parameters=paraa,description=describe)
    
    #constructing flowframe with only the non margin events
    exx2 <- exx[exx[,'Margin.Indicator']==1,]
    rflowframe <- new("flowFrame",exprs=exx2,parameters=paraa,description=describe)
  } else if(type=="estimate"){
    bfilter <- density(exprs(flow.frame)[,Channel])
    infl <- c(FALSE, diff(diff(bfilter$y)>0)>0)
    if(length(unique(infl)) > 1){
      infl_point <- max(bfilter$x[infl],na.rm=T)
    }else{
      infl_point <- max(bfilter$x,na.rm=T)
    }
    
    margin.ind <- ifelse(exprs(flow.frame)[,Channel] <= infl_point,T,F) #FALSE=Margin Event
    n_margin <- sum(margin.ind == F) #part of output, i.e. percentage
    #that are not margin events
    
    #constructing full flow frame with both margin and non-margin events
    exx <- as.matrix(cbind(exprs(flow.frame),as.numeric(margin.ind))) #1=not margin, 0=margin
    colnames(exx)[ncol(exx)] <- 'Margin.Indicator' #1=not amrgin, 0=margin
    
    #constructing the annotated data frame for the parameter
    ddata <- data.frame(rbind(flow.frame@parameters@data,
                              c("Margin.Indicator","Margin Indicator",1,0,1)))
    row.names(ddata) <- c(row.names(flow.frame@parameters@data),paste0("$","P",ncol(exx)))
    dvarMetadata <- flow.frame@parameters@varMetadata
    ddimnames <- flow.frame@parameters@dimLabels
    paraa <- Biobase::AnnotatedDataFrame(data=ddata,varMetadata=dvarMetadata,dimLabels=ddimnames)
    
    describe <- flow.frame@description
    fflowframe <- new("flowFrame",exprs=exx,parameters=paraa,description=describe)
    
    #constructing flowframe with only the non margin events
    exx2 <- exx[exx[,'Margin.Indicator']==1,]
    rflowframe <- new("flowFrame",exprs=exx2,parameters=paraa,description=describe)
  } else stop("Error: check your inputs")
  
  return(list(reducedflowframe=rflowframe,fullflowframe=fflowframe,N_margin=n_margin,
              N_cell=nrow(flow.frame)))
}
