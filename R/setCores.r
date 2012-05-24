#
#   This function is a slightly modified version of
#   function setCores in package spMC version 0.2.2
#   written by Luca Sartore <drwolf85@gmail.com>
#
setCores <-
function(n) {
  # Set number of CPU cores which will be used by the package
  #
  #     n number of CPU cores

  if (!missing(n)) {
    if (is.numeric(n)) {
      n <- as.integer(ceiling(n))
      n <- .C('setNumThreads', n = as.integer(n), DUP = FALSE, PACKAGE = "aws")$n
    }
  }
  n <- 0L
  crTot <- 0L
  n <- .C('getNumThreads', n = as.integer(n), DUP = FALSE, PACKAGE = "aws")$n
  if (n <= 1L) {
    cat("Parallel computation will not perform. CPU cores in use: 1.\n")
  }
  else {
    crTot <- .C('getNumCores', n = as.integer(crTot), DUP = FALSE, PACKAGE = "aws")$n
    cat("Parallel computation will perform.\n")
    cat("  Total CPU cores available: ", crTot, ".\n", sep = "")
    cat("  CPU cores in use: ", n, ".\n", sep = "")
  }
  invisible(n)
}
