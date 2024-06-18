#' ZLS Test for Variational Additive Models
#'
#' Performs a global test for smoothed or functional effects from variational additive models using a test for semiparametric effects.
#'
#' @param model an object of class `vam` which is returned by the function [vam()].
#' @param test an [integer] value to extract a specific test when multiple smoothed and
#'  functional effects are included in the model.
#'
#' @details Zhang, Lin, and Sowers (2000) describe a test for semiparametric effects. Meyer and Wei (2024) show that CAVI
#'  algorithm admit a form that can be used to implement this test on smoothed and functional effects estimated using
#'  mean field variational inference. For smoothed effects, the test statistic has the form
#'  \deqn{G(\textbf{Y}) = \int_{z_{(1)}}^{z_{(N)}} \textbf{Y}'\textbf{c}(z)\textbf{c}(z)'\textbf{Y} dz = \textbf{Y}'U \textbf{Y},}
#'  where \eqn{U = \int_{z_{(1)}}^{z_{(N)}} \textbf{c}(z)\textbf{c}(z)' dz}, \eqn{z_{(1)} = \min_i(z_i)}, and \eqn{z_{(N)} = \max_i(z_i)}.
#'  For functional effects, the test statistic is
#'  \deqn{G(\textbf{Y}) = \int_{t_1}^{t_T} \textbf{Y}'\textbf{c}(t)\textbf{c}(t)'\textbf{Y} dt = \textbf{Y}'Q \textbf{Y},}
#'  where \eqn{Q = \int_{t_1}^{t_T} \textbf{c}(t)\textbf{c}(t)' dt}. In either case, \eqn{G(\textbf{Y})} follows a chi-squared distribution.
#'
#' k = as.numeric(k), psi = as.numeric(psi), e = as.numeric(e),
#' comp = comp)
#' @return A list with class "`htest`" containing the following components:\tabular{ll}{
#'    \code{statistic} \tab the value of the test statistic with a name describing it. \cr
#'    \tab \cr
#'    \code{parameters} \tab  the parameters for the distribution of the test statistic. \cr
#'    \tab \cr
#'    \code{p.value} \tab the p-value for the test. \cr
#'    \tab \cr
#'    \code{null.value} \tab the location parameter. \cr
#'    \tab \cr
#'    \code{alternative} \tab a character string describing the alternative hypothesis. \cr
#'    \tab \cr
#'    \code{data.name} \tab a character string giving the names of the data. \cr
#'    \tab \cr
#'    \code{method} \tab character string giving the type of doubly ranked test performed. \cr
#'    \tab \cr
#'    \code{k} \tab scaling term for test statistic. \cr
#'    \tab \cr
#'    \code{psi} \tab approximately equal to \eqn{2tr\{ (QV)^2 \}} or \eqn{2tr\{ (UV)^2 \}}, used for degrees of freedom calculation. \cr
#'    \tab \cr
#'    \code{e} \tab approximately equal to \eqn{tr(QV) } or \eqn{tr(UV) }, used for degrees of freedom calculation. \cr
#' }
#'
#' @references Zhang, D, Lin, X, and Sowers, M (2000). Semiparametric regression for periodic longitudinal hormone data from multiple menstrual cycles.
#'  \emph{Biometrics} **56**: 31–39.
#'
#'  Zhang, D, and Lin, X (2003). Hypothesis testing in semiparametric additive mixed models. \emph{Biostatistics} **4**: 57–74.
#'
#'  Meyer, MJ and Wei, J (2024). Global tests for smoothed functions in mean field variational additive models. \emph{Available on arXiv} at \url{https://arxiv.org/abs/2406.08168}
#'
#' @examples
#' ##### Semiparametric Effect #####
#' library(SemiPar)
#' data(lidar)
#'
#' modLi   <- vam(fixed = logratio ~ 1, se = ~ range, data = lidar,
#'                se.spline = list(df = 8, degree = 3, slength = 100))
#' zls_test(modLi)
#'
#' ##### Functional Effect #####
#' library(refund)
#' DTI2_1v         <- DTI2[which(DTI2$visit == 1), ]
#' DTI2_1v$pasatb  <- 1*(DTI2_1v$pasat > median(DTI2_1v$pasat))
#' cca             <- refund::fpca.sc(DTI2_1v$cca)
#' fcca            <- cca$Yhat
#'
#' ###### CCA ######
#' modCc   <- vam(fixed = pasatb ~ 1, fe = ~ fcca, data = DTI2_1v, Ts = list(fcca = ncol(fcca)),
#'                family = binomial(link = "probit"),
#'                fe.spline = list(df = 9, degree = 3, xi = 0.01))
#' zls_test(modCc)
#'
#' @export
zls_test <- function(model, test = NULL){
  UseMethod("zls_test")
}
#' @rdname zls_test
#' @export
zls_test.vam <- function(model, test = NULL){

  family      <- model$model$family$family
  Th          <- model$model$Th
  SigmaMat    <- model$varComps$Sigmaq
  p           <- model$model$p
  Kf          <- model$model$Kf
  Kfi         <- model$model$Kfi
  Ctf         <- model$Ctf
  BqeVec      <- model$varComps$Bqe
  Muqb        <- model$VBEstimates[c(1:model$model$p)]
  Muqa        <- model$Muqa
  Ae          <- model$priors$Ae
  Y           <- model$model$Y
  test.flag   <- model$model$test.flag
  n           <- length(Y)
  BqfVec      <- model$varComps$Bqf
  XZ          <- model$model$dmat$Xc
  out         <- vector('list', length = length(Th))
  for(f in 1:length(Th)){

    ## ZLS test ##
    Sigmaqf   <- SigmaMat[(Kfi[f] + 1):(sum(Kfi[f:(f+1)])),(Kfi[f] + 1):(sum(Kfi[f:(f+1)]))]
    ctf       <- Th[[f]]%*%Ctf[(Kfi[f] + 1):(sum(Kfi[f:(f+1)])),]

    if(family == "gaussian"){
      V0        <- ((BqeVec[1] + sum(BqfVec))/(Ae + n/2 - 1))*diag(n)
    } else{
      XZMuqa    <- XZ%*%Muqb
      Vau0      <- 1 + (XZMuqa*stats::dnorm(-XZMuqa)/stats::pnorm(-XZMuqa)) - (stats::dnorm(-XZMuqa)/stats::pnorm(-XZMuqa))^2
      Vau1      <- 1 - (XZMuqa*stats::dnorm(-XZMuqa)/(1 - stats::pnorm(-XZMuqa))) - (stats::dnorm(-XZMuqa)/(1 - stats::pnorm(-XZMuqa)))^2
      Vars      <- ((1*(sign(Muqa) < 0))*Vau0 + (1*(sign(Muqa) >= 0))*Vau1)
      V0        <- diag(as.numeric(Vars))
    }
    Ct        <- t(ctf)%*%ctf

    S         <- t(Y)%*%Ct%*%Y
    CV        <- Ct%*%V0
    e         <- sum(diag(CV))
    psi       <- 2*(sum(diag(CV%*%CV)))
    k         <- psi/(2*e)
    nu        <- 2*(e^2)/psi
    Xobs      <- S/k
    pval      <- 1-stats::pchisq(Xobs, df = nu)

    X2              <- as.numeric(Xobs)
    names(X2)       <- 'Chi-squared'
    nullval         <- 0
    if(test.flag[f] == 'se'){
      names(nullval)  <- 's(x)'
    } else{
      names(nullval)  <- 'g(t)'
    }

    comp      <- list(S = S, Ctf = Ctf, Ct = Ct, V0 = V0, ct = ctf,
                      Sigmaqf = Sigmaqf)
    out[[f]]  <- list(statistic = X2, p.value = as.numeric(pval),
                      parameters = list(df = as.numeric(nu)),
                      method = 'One Sample Variational ZLS Test',
                      alternative = 'two-sided', null.value = nullval,
                      data.name = model$data.name[[f]],
                      k = as.numeric(k), psi = as.numeric(psi), e = as.numeric(e),
                      comp = comp)
    class(out[[f]])   <- 'htest'
  }

  if(length(out) > 1){
    if(is.null(test)){
      return(out)
    } else{
      return(out[[test]])
    }
  } else{
    return(out[[1]])
  }
}
