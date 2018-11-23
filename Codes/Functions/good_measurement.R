# indicates if measurement from a flowfile is good or bad
#
# @param fcsfile flowfile to be checked.
# @param metafile associated metafile to the supplied fcsfile.
# @param CpermL column name in metafile containing cell per microlitre.
# @param rowinmetafile row containing measurements associated to fcsfile.
# @param d_cellnum desired minimum number of particle counts. Flowfiles with smaller particle counts are
#                  termed bad. Defaults to 90.
# @param d_cellpML maximal accepted cell per microlitre. Flowfiles with larger cell per microlitre are
#        termed bad. Defaults to 1000
# @return dataframe with columns
# \itemize{
#         \item \strong{CML -} cell per microlitre as read from the associated metafile.
#         \item \strong{status -} flowfile status, either "good" or "bad".
#         \item \strong{ID -} is the row number in meta file
#         \item \strong{Size -} is the number of cells in the FCS file.
#         }
# @examples
# \dontrun{
#   goodfcs(fcsfile = flowframe, metadtat, Cell_per_microlitre, 12, d_cellnum = 500, 1500)
# }
#
# @export goodfcs

goodfcs <- function(fcsfile, metafile, CpermL, rowinmetafile, d_cellnum = 90, d_cellpML = 1000){
  if(!is.null(metafile) & !is.null(fcsfile) & !is.null(CpermL) & !is.null(rowinmetafile)){
    cML <- as.numeric(metafile[rowinmetafile, CpermL])
    size <- nrow(fcsfile)
    goodfile <- ifelse((cML < d_cellpML) & (size > d_cellnum),"good","bad")
  } else stop("One or more of the inputs is empty")
  return(data.frame(CML = cML, status = goodfile, ID = rowinmetafile, Size = size))
}
