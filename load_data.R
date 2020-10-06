#' Import protein or metabolite dataset
#'
#' This function loads your dataset containing elution profiles of proteins or metabolites.
#'
#' @param x Path to the input file
#' @return A data frame with elution profiles
#' @export

import_data <- function(x, nr_replicas){
  #meta_data <- read.delim(paste("H:/Software, devices, machines/R/PROMIS package/db/", x, sep = ""), row.names = 1)
  meta_data <- read.delim(x)
  i <- 1
  while(i <= nr_replicas){
    tmp1 <- length(grep(paste0("Rep_", i), colnames(meta_data)))
    print(paste("Number of fractions in Rep_", i, " = ", tmp1, sep = ""))
    i <- i+1
  }

  meta_data
}
