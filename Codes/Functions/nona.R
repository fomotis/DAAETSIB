# Removes NA values from the expression matrix of a flowfile
#
# @param x is the flowframe with expression matrix containing NA
# @return a flowframe with expression matrix rid of NAs
#

nona <- function(x){
  dtest <- !apply(flowCore::exprs(x), 1, function (row) any(is.na(row) | is.nan(row)))
  exx <- flowCore::exprs(x)[dtest == T, ]
  paraa <- x@parameters
  describe <- x@description
  new("flowFrame", exprs = exx, parameters = paraa, description = describe)
}
