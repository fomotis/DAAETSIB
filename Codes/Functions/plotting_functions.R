pair_plot <- function(x,notToPlot=c("TIME")){
  col.palette <- colorRampPalette(c("white","blue","cyan","green","orange", "red"),
                                  space = "Lab")
  pairs(exprs(x)[,which(!(colnames(x) %in% notToPlot ))],pch=".",
        panel = function(...) smoothScatter(..., nrpoints = 0, colramp = col.palette,
                                            add = TRUE),gap=0.2,
        main=identifier(x))
}