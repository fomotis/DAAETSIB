# Removes negative values from the expression matrix
#
#@param x is the flowframe whose expression matrix contains negative values
#@return a flow frame with non-negative values in its expression matrix
#
noneg <- function(x){
  dtest <- !apply(flowCore::exprs(x), 1 ,function (row) any(row <= 0))
  exx <- flowCore::exprs(x)[dtest == T, ]
  paraa <- x@parameters
  describe <- x@description
  methods::new("flowFrame", exprs = exx, parameters = paraa, description = describe)
}
