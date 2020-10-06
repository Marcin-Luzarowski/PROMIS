#' Plot elution profile of a molecule
#' @param x Data frame containing single elution profiles
#' @param met_plot_color Color of a plot
#' @param do_plot If plot should be made
#' @param normalize If data should be maximum transformed
#' @param tmp_dir Determine temporary directory to which the peaks will be saved
#' @return A data frame without elution profiles, which peaks span less than c fractions
#' @export

plot_single_profiles <- function(x, met_plot_color, do_plot, normalize, tmp_dir, fraction_names){
  rownames(x) <- x$Name
  x <- x[,grep("Fraction", colnames(x))]
  if(normalize == TRUE){
    x <- maxNormalize(x)
  }
  if(do_plot == TRUE){
    print(paste("Elution profiles are now stored in ", tmp_dir, sep = ""))
    i <- 1
    while(i <= nrow(x)){
      plot_name <- paste(tmp_dir, "/",  rownames(x)[i], ".jpeg", sep = "")
      jpeg(plot_name, width = 8, height = 8, units = 'cm', res = 600)
      par(mfrow = c(1,1), xpd = TRUE, mar=c(4,2.5,3,1.5))
      plot(1:nr_fractions, x[i,], col=met_plot_color, type = "l", ylim = c(0,1), lwd = 1, pch = 0, axes = FALSE, ann=FALSE, lty = 1)
      legend("topleft", inset = c(-0, -0.2), rownames(x)[i], col=met_plot_color, lty = 1, cex = 0.5)
      axis(1, at=1:nr_fractions, lab=fraction_names, las = 2, cex.axis = 0.3)
      axis(2, las=1, cex.axis = 0.4)
      box()
      title(ylab = list("Relative intensity (Max)", cex = 0.5), line = 1.5)
      title(xlab = list("Fraction number", cex = 0.5))
      dev.off()
      i <- i+1
    }
  }
  return(x)
}
