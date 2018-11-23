### Function to determine if the current FCS file is good for further use.
#CpermL is the column containing cell per micro litre in metafile from flow cytometer
#rowinmetafile is the row containing information corresponding to the fcsfile
#d_cellnum is least the desired number of cells needed
#d_cellpML is the desired cell per micor litre


goodfcs <- function(fcsfile,metafile,CpermL,rowinmetafile,d_cellnum=90,d_cellpML=1000){
  if(!is.null(metafile) & !is.null(fcsfile) & !is.null(CpermL) & !is.null(rowinmetafile)){
    cML <- as.numeric(metafile[rowinmetafile,CpermL])
    size <- nrow(fcsfile)
    goodfile <- ifelse((cML < d_cellpML) & (size > d_cellnum),"good","bad")
  } else stop("One or more of the inputs is empty")
  return(data.frame(CML=cML,status=goodfile,ID=rowinmetafile,Size=size))
}