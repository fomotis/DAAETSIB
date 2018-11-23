# log transforms the expression matrix of a flowframe
#
# @param x flowframe to be transformed
# @param notToTransform columns not to be transformed
#
# @return \strong{flowframe} with log transformed expression matrix
#
# @examples
#
# \dontrun{
#lnTrans(margin_filtered$Day1[[1]], c("YEL.B.HLin", "TIME"))
# }
#
# @export lnTrans

lnTrans <- function(x, notToTransform = c("SSC.W", "TIME")){
  exx <- cbind(log(flowCore::exprs(x)[, which(!(colnames(x) %in% notToTransform ))]),
               flowCore::exprs(x)[, which(colnames(x) %in% notToTransform)])
  colnames(exx) <- colnames(x)
  paraa <- x@parameters
  describe <- x@description
  methods::new("flowFrame", exprs = exx, parameters = paraa, description = describe)
}
