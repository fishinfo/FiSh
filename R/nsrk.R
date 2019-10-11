#' Normal scale rule for kernel density estimation
#'
#' Bandwidth selector for non-parametric estimation. Estimates
#' the optimal AMISE bandwidth using the Normal Scale Rule with Gaussian kernel.
#'
#' @usage nsrk(x, log_trsf=FALSE)
#' @param x  Univariate data.
#' @param log_trsf Logical flag: if \code{TRUE} the data are log-transformed (usually used for skewed positive data).
#' By default \code{log_trsf = FALSE}.
#'
#' @return The bandwidth value.
#'
#'
#' @examples
#'
#' x <- rnorm(1000)
#' h <- nsrk(x)
#'
#' @references
#' M. P. Wand and M. C. Jones,  (1995). Kernel Smoothing. Chapman and Hall, London.
#'
#' @importFrom stats quantile var
#' @export
#'
# Normal scaled rule
nsrk <- function(x, log_trsf=FALSE){
  n <- length(x)
  if (log_trsf){
    if (any(x<0)) {stop("Data must be positive when log_trsf is TRUE.")}
    x <- log(x)
  }
  scalest <- min(stdev = sqrt(var(x)), iqr = (quantile(x,3/4) - quantile(x, 1/4))/1.349)

  return(scalest*((4/(3*n))^(1/5)))
}
