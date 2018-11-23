### Output
# fullframe = the flowframe supplied (which should be the reduced frame from the method above)
# rframe = the flowframe after the doublets have been filterd out
# proportion is the percentage of singlets with respect to rframe
# indicator is an indicator as to which cells (rows) are singlets and which are doublets
#  1 = singlet, 0=doublet
#### Fishing out Doublets
celldoublets <- function(rframe,p1="SSC.HLin",p2="SSC.ALin",div=NULL,
                         manualgatetype=c("rectangle","polygon","ellipsoid"),
                         type=c("manual","flowDens"),positions=c(F,F),control=NULL,
                         control_use=c(T,T),control_flowFrame=NULL){
  if(type=="manual"){
    #column names of div
    if(manualgatetype=="rectangle" & !is.null(div)){
      colnames(div) <- c(p1,p2)
      drg <- rectangleGate(filterId="rectid",.gate=div)
      drest <- flowCore::filter(rframe,drg)
      dind <- drest@subSet
      n_singlet <- sum(dind)
      n_particles <- nrow(rframe) #part of output
      #flowframe with doublets filtered out
      dreducedframe <- flowCore::Subset(rframe,drg) #part of output
      #plotting
      plot(exprs(rframe)[,p1],exprs(rframe)[,p2],main="Doublets Gating",
           xlab=p1,ylab=p2,pch=".",col=c("red3","green3")[factor(dind)],cex=1.5,frame.plot=F)
      legend("bottomright",c("singlets","doublets"),col=c("green3","red3"),pch=c(".","."),
             cex=1.5,bty="n")
      text(4,6,labels=paste0("ID=",counter))
      
      #flowframe with indicator for doublets
      exx <- cbind(exprs(rframe),as.numeric(dind)) #1=singlet, 0=doublet
      colnames(exx)[length(colnames(exx))] <- 'Singlet.Indicator' #1=singlet, 0=doublet
      
      #constructing the annotated data frame for the parameter
      rframe@parameters@data[,1] <-
        as.character(rframe@parameters@data[,1])
      ddata <- data.frame(rbind(rframe@parameters@data,
                                c("Singlet.Indicator","Singlet Indicator",1,0,1)))
      row.names(ddata) <- c(row.names(rframe@parameters@data),'$P14')
      dvarMetadata <- rframe@parameters@varMetadata
      ddimnames <- rframe@parameters@dimLabels
      paraa <- Biobase::AnnotatedDataFrame(data=ddata,varMetadata=dvarMetadata,
                                           dimLabels=ddimnames)
      describe <- rframe@description
      fflowframe <- new("flowFrame",exprs=exx,parameters=paraa,description=describe)
    }else if(manualgatetype=="polygon" & !is.null(div)){
      colnames(div) <- c(p1,p2)
      dpg <- polygonGate(filterid="dpg",.gate=div)
      drest <- flowCore::filter(rframe,dpg)
      dind <- drest@subSet
      n_singlet <- sum(dind)
      n_particles <- nrow(rframe)
      #flowframe with doublets filtered out
      dreducedframe <- flowCore::Subset(rframe,drest) #part of output
      
      #plotting
      plot(exprs(rframe)[,p1],exprs(rframe)[,p2],main="Doublets Gating",
           xlab=p1,ylab=p2,pch=".",col=c("red3","green3")[factor(dind)],cex=1.5,frame.plot=F)
      legend("bottomright",c("singlets","doublets"),col=c("green3","red3"),pch=c(".","."),
             cex=1.5,bty="n")
      text(4,6,labels=paste0("ID=",counter))
      
      #flowframe with indicator for doublets
      exx <- cbind(exprs(rframe),as.numeric(dind)) #1=singlet, 0=doublet
      colnames(exx)[length(colnames(exx))] <- 'Singlet.Indicator' #1=singlet, 0=doublet
      
      #constructing the annotated data frame for the parameter
      rframe@parameters@data[,1] <-
        as.character(rframe@parameters@data[,1])
      ddata <- data.frame(rbind(rframe@parameters@data,
                                c("Singlet.Indicator","Singlet Indicator",1,0,1)))
      row.names(ddata) <- c(row.names(rframe@parameters@data),'$P14')
      dvarMetadata <- rframe@parameters@varMetadata
      ddimnames <- rframe@parameters@dimLabels
      paraa <- Biobase::AnnotatedDataFrame(data=ddata,varMetadata=dvarMetadata,
                                           dimLabels=ddimnames)
      describe <- rframe@description
      fflowframe <- new("flowFrame",exprs=exx,parameters=paraa,description=describe)
    }else if(manualgatetype=="ellipsoid"){
      ms <- c(mean(exprs(rframe)[,p1]),mean(exprs(rframe)[,p2]))
      cova <- cov(exprs(rframe)[,c(p1,p2)])
      deg <- ellipsoidGate("egid",.gate=cova,mean=ms,distance=4)
      drest <- flowCore::filter(rframe,deg)
      dind <- drest@subSet
      n_singlet <- sum(dind)
      n_particles <- nrow(rframe)
      #flowframe with doublets filtered out
      dreducedframe <- flowCore::Subset(rframe,drest) #part of output
      
      #plotting
      plot(exprs(rframe)[,p1],exprs(rframe)[,p2],main="Doublets Gating",
           xlab=p1,ylab=p2,pch=".",col=c("red3","green3")[factor(dind)],cex=1.5,frame.plot=F)
      legend("bottomright",c("singlets","doublets"),col=c("green3","red3"),pch=c(".","."),
             cex=1.5,bty="n")
      text(4,6,labels=paste0("ID=",counter))
      
      #flowframe with indicator for doublets
      exx <- cbind(exprs(rframe),as.numeric(dind)) #1=singlet, 0=doublet
      colnames(exx)[length(colnames(exx))] <- 'Singlet.Indicator' #1=singlet, 0=doublet
      
      #constructing the annotated data frame for the parameter
      rframe@parameters@data[,1] <-
        as.character(rframe@parameters@data[,1])
      ddata <- data.frame(rbind(rframe@parameters@data,
                                c("Singlet.Indicator","Singlet Indicator",1,0,1)))
      row.names(ddata) <- c(row.names(rframe@parameters@data),'$P14')
      dvarMetadata <- rframe@parameters@varMetadata
      ddimnames <- rframe@parameters@dimLabels
      paraa <- Biobase::AnnotatedDataFrame(data=ddata,varMetadata=dvarMetadata,
                                           dimLabels=ddimnames)
      describe <- rframe@description
      fflowframe <- new("flowFrame",exprs=exx,parameters=paraa,description=describe)
    }else stop("Error: Select a manual gate type and provide the corresponding arguments")
    
  }else if(type=="flowDens" & control==F){
    tfden2 <- flowDensity(rframe,channels=c(p1,p2),position=positions,ellip.gate=T,
                          upper=c(T,T))
    plot(getflowFrame(tfden2),tfden2)
    text(tfden2@gates[1],tfden2@gates[2],labels=paste0("ID=",counter))
    
    
    ttt1 <- exprs(tfden2@flow.frame) #with NA values for doublets
    dind <- ifelse(!is.na(ttt1[,1])==T,T,F) #!NA, i.e. T is a singlet
    n_singlet <- sum(dind)
    n_particles <- nrow(rframe)
    exx <- cbind(exprs(rframe),as.numeric(dind)) #1=singlet, 0=doublet
    colnames(exx)[length(colnames(exx))] <- 'Singlet.Indicator' #1=singlet, 0=doublet
    
    #constructing the annotated data frame for the parameter
    rframe@parameters@data[,1] <-
      as.character(rframe@parameters@data[,1])
    ddata <- data.frame(rbind(rframe@parameters@data,
                              c("Singlet.Indicator","Singlet Indicator",1,0,1)))
    row.names(ddata) <- c(row.names(rframe@parameters@data),'$P14')
    dvarMetadata <- rframe@parameters@varMetadata
    ddimnames <- rframe@parameters@dimLabels
    paraa <- Biobase::AnnotatedDataFrame(data=ddata,varMetadata=dvarMetadata,
                                         dimLabels=ddimnames)
    
    describe <- rframe@description
    fflowframe <- new("flowFrame",exprs=exx,parameters=paraa,description=describe)
    
    #flowframe with doublets filtered out
    dreducedframe <- getflowFrame(tfden2) # part of the output
  } else if(type=="flowDens" & control==T & !is.null(control_flowFrame)){
    tfden2 <- flowDensity(rframe,channels=c(p1,p2),position=positions,use.control=control_use,
                          control=control_flowFrame,ellip.gate=T,upper=c(T,T))
    plot(getflowFrame(tfden2),tfden2)
    text(tfden2@gates[1],tfden2@gates[2],labels=paste0("ID=",counter))
    autoplot(getflowFrame(tfden2),p1,p2)
    ttt1 <- exprs(tfden2@flow.frame) #with NA values for doublets
    dind <- ifelse(!is.na(ttt1[,1])==T,T,F) #!NA, i.e. T is a singlet
    n_singlet <- sum(dind)
    n_particles <- nrow(rframe)
    exx <- cbind(exprs(rframe),as.numeric(dind)) #1=singlet, 0=doublet
    colnames(exx)[length(colnames(exx))] <- 'Singlet.Indicator' #1=singlet, 0=doublet
    
    #constructing the annotated data frame for the parameter
    rframe@parameters@data[,1] <-
      as.character(rframe@parameters@data[,1])
    ddata <- data.frame(rbind(rframe@parameters@data,
                              c("Singlet.Indicator","Singlet Indicator",1,0,1)))
    row.names(ddata) <- c(row.names(rframe@parameters@data),'$P14')
    dvarMetadata <- rframe@parameters@varMetadata
    ddimnames <- rframe@parameters@dimLabels
    paraa <- Biobase::AnnotatedDataFrame(data=ddata,varMetadata=dvarMetadata,
                                         dimLabels=ddimnames)
    
    describe <- rframe@description
    fflowframe <- new("flowFrame",exprs=exx,parameters=paraa,description=describe)
    dprop <- tfden2@proportion
    #flowframe with doublets filtered out
    dreducedframe <- getflowFrame(tfden2) # part of the output
  }
  else stop("Select a GateType")
  return(list(fullframe=fflowframe,reducedframe=dreducedframe,N_Singlets=n_singlet,
              N_Cells=n_particles, singletindicator=as.numeric(dind)))
}