
debris_nc <- function(flowframe, p1, p2) {
  #debris gating
  if(diptest::dip.test(flowframe@exprs[, p1])$p.value < 0.05){
    gate_deb <- flowDensity::deGate(flowframe, p1, bimodal = T)
  } else gate_deb <- flowDensity::deGate(flowframe, p1, bimodal = T, upper = F)
  #plotting
  graphics::par(mfrow = c(1, 2))
  flowDensity::plotDens(flowframe, c(p1, p2), xlab = p1, ylab = p2,
                        main = paste0("ID = ", flowCore::identifier(flowframe)), 
                        frame.plot = F)
  abline(v = gate_deb, lty = 3, lwd = 2, col = 2)
  debs_exp <- flowCore::exprs(flowframe)[flowCore::exprs(flowframe)[, p1] <= gate_deb, ]
  text(mean(debs_exp[, p1]), mean(debs_exp[, p2]), "Deb", col = 2)
  
  #reduced_flowframe containing BS4 and BS5 only
  bs4bs5 <- flowframe[which(flowCore::exprs(flowframe)[, p1] > gate_deb), ]
  #positions
  deb_pos <- which(flowCore::exprs(flowframe)[, p1] <= gate_deb)
  other_pos <- which(flowCore::exprs(flowframe)[, p1] > gate_deb)
  return(list(bs4bs5 = bs4bs5, deb_pos = deb_pos, other_pos = other_pos))
}
