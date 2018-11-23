# function to plot the expression matrix of a flowframe. Note that, it takes some time to
# display the plot.
#
# @param x flowframe to be plotted
# @param notToPlot column in expression matrix not to be plotted
#
# @examples
# \dontrun{
#   pair_plot(log_transformedSet$Day1[[1]])
# }
#
#
# @importFrom grDevices colorRampPalette
# @importFrom graphics abline points panel.smooth pairs smoothScatter text
# @export pair_plot
#

pair_plot <- function(x, notToPlot = c("TIME")) {
  col.palette <- colorRampPalette(c("white","blue","cyan","green","orange", "red"),
                                  space = "Lab")
  pairs(flowCore::exprs(x)[, which(!(colnames(x) %in% notToPlot ))], pch = ".",
        panel = function(...) smoothScatter(..., nrpoints = 0, colramp = col.palette,
                                            add = TRUE), gap = 0.2,
        main = flowCore::identifier(x))
}
