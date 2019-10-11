#' Fisher-Shannon method
#'
#' Non-parametric estimates of the Shannon Entropy Power (SEP), the Fisher Information Measure (FIM) and the
#' Fisher-Shannon Complexity (FSC), using kernel density estimators with Gaussian kernel.
#'
#' @usage SEP_FIM(x, h, log_trsf=FALSE, resol=1000, tol = .Machine$double.eps)
#' @param x  Univariate data.
#' @param h Value of the bandwidth for the density estimate
#' @param log_trsf Logical flag: if \code{TRUE} the data are log-transformed (used for skewed data), in this case
#' the data should be positive. By default, \code{log_trsf = FALSE}.
#' @param resol Number of equally-spaced points, over which function approximations are computed and integrated.
#' @param tol A tolerance to avoid dividing by zero values.
#'
#' @return A table with one row containing:
#'  \itemize{
#'   \item \code{SEP} Shannon Entropy Power.
#'   \item \code{FIM} Fisher Information Measure.
#'   \item \code{FSC} Fisher-Shannon Complexity
#'   }
#'
#' @examples
#'
#' library(KernSmooth)
#' x <- rnorm(1000)
#' h <- dpik(x)
#' SEP_FIM(x, h)
#'
#'
#'
#' @references
#' F. Guignard, M. Laib, F. Amato, M. Kanevski, Advanced analysis of temporal
#' data using Fisher-Shannon information : theoretical development and
#' application to geoscience
#'
#' @import fda.usc KernSmooth
#' @export

SEP_FIM <- function(x, h, log_trsf=FALSE, resol=1000, tol = .Machine$double.eps){
  if (log_trsf){
    if(any(x<0)){stop("Data must be positive when log_trsf is TRUE.")}
    return(SEP_FIM_log(x, h, tol=tol))
  } else{
    return(SEP_FIM_nolog(x, h, tol=tol))
  }

}

SEP_FIM_nolog <- function(x,  h, resol=1000, tol = .Machine$double.eps){
  n <- length(x)
  integ_start<-min(x)
  integ_end<-max(x)
  Xgrid<-seq(integ_start, integ_end, length.out = resol)
  Accu_f<-rep(0, resol)
  Accu_f_drv<-rep(0, resol)
  for (i in 1:n){
    ddist<-Xgrid-x[i]
    kern <- exp((-ddist^2)/(2*h^2))
    Accu_f <- Accu_f + kern
    Accu_f_drv <- Accu_f_drv + ddist*kern
  }


  FIM<- ifelse(Accu_f > tol, Accu_f_drv^2/Accu_f, 0)

  FIM <- int.simpson(fdata(FIM,Xgrid))
  FIM <- FIM*(1/(sqrt(2*pi)*n*h^5))
  # SEP
  f <- Accu_f/(sqrt(2*pi)*n*h)
  H <- ifelse(f>tol, -f*log(f), 0)
  H <- int.simpson(fdata(H,Xgrid))
  SEP<-(1/(2*pi*exp(1)))*exp(2*H)

  return(data.frame(SEP=round(SEP,4), FIM=round(FIM,4), FSC=round(SEP*FIM,4)))
}


SEP_FIM_log<-function(x,  h, resol=1000, tol = .Machine$double.eps){
  n <- length(x)
  integ_start<-min(x)
  integ_end<-max(x)
  # Transform :
  y <- log(x)
  Ygrid <- seq(log(integ_start), log(integ_end), length.out = resol)
  Accu_f<-rep(0, resol)
  Accu_f_drv<-rep(0, resol)
  for (i in 1:n){
    ddist<-Ygrid-y[i]
    kern <- exp(-(ddist^2)/(2*h^2))
    Accu_f <- Accu_f + kern
    Accu_f_drv <- Accu_f_drv + ddist*kern
  }
  back_Xgrid <- exp(Ygrid)
  Accu_f <- Accu_f/h
  Accu_f_drv <- (Accu_f_drv/h^3)+(Accu_f)

  #FIM <- (Accu_f_drv^2/Accu_f)
  FIM<- ifelse(Accu_f > tol, Accu_f_drv^2/Accu_f, 0)
  FIM <- FIM/back_Xgrid^3
  #FIM<-integrate.xy(back_Xgrid, FIM)
  FIM <- int.simpson(fdata(FIM,back_Xgrid), equi = F)
  FIM <- FIM*(1/(sqrt(2*pi)*n))

  f <- Accu_f/(sqrt(2*pi)*n*back_Xgrid)
  #H <- -f*log(f)
  #H<-integrate.xy(back_Xgrid, H)
  H <- ifelse(f>tol, -f*log(f), 0)
  H <- int.simpson(fdata(H,back_Xgrid), equi = F)
  SEP<-(1/(2*pi*exp(1)))*exp(2*H)

  return(data.frame(SEP=round(SEP,4), FIM=round(FIM,4), FSC=round(SEP*FIM,4)))
}
