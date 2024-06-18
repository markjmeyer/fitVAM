#' Variational Additive Model
#'
#' Implements variational algorithms for estimating effects from both Gaussian and binary additive models with
#' arbitrary numbers of smoothed/semiparametric and functional effects.
#'
#' @param fixed an object of class [formula] containing the fixed effects.
#' @param se a one-sided object of class [formula] containing the smoothed effects.
#' @param fe a one-sided object of class [formula] containing the functional effects.
#' @param Ts a vector containing the grid on which the function is measured.
#' @param data [data.frame] containing the variables in the fixed effects model.
#' @param se.spline a list containing three objects specifying the spline for the smoothed effects: `df` degrees of freedom,
#'  `degree` of the spline, and `slength` the number of evaluation points to use when graphing the smoothed effect.
#'  Defaults to `list(df = 8, degree = 3, slength = 100)`.
#' @param fe.spline a [list] containing three objects specifying the spline for the functional effect:  `df` degrees of freedom,
#'  `degree` of the spline, and `xi` a sparsity control between 0 and 1. Defaults to `list(df = 16, degree = 3, xi = 0.01)`.
#' @param family a description of the error distribution and link function to be used in the model.
#' Currently only `gaussian(link = "identity")` and `binomial(link = "probit")` are supported.
#' @param priors a [list] containing the prior specifications for the model errors (`Ae` and `Be`),
#'  the functional or smoothed effects (`Af` or `Bf`, can be vectors to vary priors on smoothed tuning parameters),
#'  and the variance of the fixed effects (`sigB`). Defaults to
#'  `list(Ae = 0.01, Be = 0.01, Af = c(0.01, 0.01), Bf = c(0.01, 0.01), sigB = 1000)`.
#' @param tol [numeric] value indicating the threshold or tolerance at which the model is judged to have converged. Defaults to 1e-5.
#' @param maxIter [numeric] value indicating the maximum number of iterations to perform should `tol` not be reached. Defaults to 500.
#' @param up an [integer] value dictating the frequency of printed iteration update statements. Defaults to 100.
#' @param dots an [integer] value dictating the frequency of printed iteration update dots. Defaults to 10.
#' @param verbose a [logical] value. If `TRUE`, iteration updates are printed in the console. If `FALSE`, they are suppressed.
#'  Defaults to `TRUE`.
#' @details Varitaional additive models are additive models fit using variational Bayesian inference. The models generally have the form
#' \deqn{g(y_i) = \alpha + \boldsymbol{\beta}'\textbf{x}_i + \sum_{m = 1}^M s_m(z_{im}) + \sum_{f = 1}^F \int_{t\in\mathcal{T}} w_{if}(t)\gamma_f(t) dt + \epsilon_i,}
#' where \eqn{s()} and \eqn{\gamma()} are smoothed and functional effects, respectively. \eqn{g()} is a link function determined by the
#' arugment `family`. Currently, `vam` only supports `family = gaussian(link = "identity")` or `family = binomial(link = "probit")`.
#' Variational inference is a deterministic, algorithmic approach to estimating Bayesian posterior densities.
#' Both smoothed and functional effects are penalized via the mixed model representation of semiparametric regression effects.
#'
#' Each model component is specified separately using the `fixed` argument for fixed effects (i.e. non-penalized effects), `se` for
#' smoothed effects, and `fe` for functional effects.
#'
#' @return A list with class "`vam`" containing the following components:\tabular{ll}{
#'    \code{VBEstimates} \tab The last iteration of variational Bayesian estimates of the mean model. \cr
#'    \tab \cr
#'    \code{Muqa} \tab  if family is binomial, the last latent vector otherise [NULL]. \cr
#'    \tab \cr
#'    \code{fixedEffects} \tab estimates of \eqn{\boldsymbol{\beta}}. \cr
#'    \tab \cr
#'    \code{smoothEffect} \tab estimates of the smoothed effects. \cr
#'    \tab \cr
#'    \code{funcEffect} \tab estimates of the functional effects. \cr
#'    \tab \cr
#'    \code{Ctf} \tab last vector \eqn{c(t)} for use in ZLS test. \cr
#'    \tab \cr
#'    \code{data.name} \tab a character string giving the names of the data. \cr
#'    \tab \cr
#'    \code{varComps} \tab a list containing model variance estimates: `Sigmaq` last covariance,
#'    `Bqe` last component of model variance, and `Bqf` last components of smoothed/functional effects. \cr
#'    \tab \cr
#'    \code{model} \tab a list containing model details including: arguments `fixed`, `fe`, `se`, `data`, and `family`
#'    as well as the number of fixed effects `p`, number of spline bases `Kf` (total) and `Kfi` (by effect), the design matrix (`dmat`),
#'    penalty matrices (`O`), basis matrices (`Th`), the outcome `Y`, and a vector indicating effect type for testing (`test.flag`). \cr
#'    \tab \cr
#'    \code{priors} \tab argument `prior`. \cr
#'    \tab \cr
#'    \code{nfun} \tab total number of functions fit, smoothed and functional. \cr
#'    \tab \cr
#'    \code{algOut} \tab list of algorithm details including: total iterations to convergence (`totalIter`), `tol`, `maxIter`, and
#'      a list (`logpd`) containing details on the changes in the K-L lower bound
#'      (`lDelta` tracks changes, `logp` tracks current lower bound, `prev_logp` tracks previous). \cr
#'    \tab \cr
#'    \code{se.spline} \tab argument `se.spline`. \cr
#'    \tab \cr
#'    \code{fe.spline} \tab argument `fe.spline`. \cr
#' }
#'
#' @references Meyer, MJ and Wei, J (2024). Global tests for smoothed functions in mean field variational additive models. \emph{Available on arXiv} at \url{https://arxiv.org/abs/2406.08168}
#'
#' @examples
#' ##### Semiparametric Effect #####
#' library(SemiPar)
#' data(lidar)
#'
#' modLi   <- vam(fixed = logratio ~ 1, se = ~ range, data = lidar,
#'                se.spline = list(df = 8, degree = 3, slength = 100))
#' plot(modLi)
#' coef(modLi)
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
#' plot(modCc)
#' coef(modCc)
#'
#' @export
vam   <- function(fixed, se = NULL, fe = NULL, Ts = NULL, data,
                  se.spline = NULL, fe.spline = NULL, family = stats::gaussian(link = "identity"),
                  priors = list(Ae = 0.01, Be = 0.01, Af = c(0.01, 0.01),
                                Bf = c(0.01, 0.01), sigB = 1000),
                  tol = 1e-5, maxIter = 500, up = 100, dots = up/10,
                  verbose = TRUE){

  if(family$family != "gaussian" & family$family != "binomial"){
    stop("vam only supports family = gaussian(link = 'identity') or family = binomial(link = 'probit')")
  }
  if(family$family == "binomial" & family$link != "probit"){
    stop("vam only supports family = binomial(link = 'probit') for binary outcomes")
  }
  if(family$family == "gaussian" & family$link != "identity"){
    stop("vam only supports family = gaussian(link = 'identity') for continuous outcomes")
  }
  ## extract data, build matrices and vectors
  mf      <- stats::model.frame(formula = fixed, data = data)
  Y       <- as.numeric(stats::model.response(mf))
  Xc      <- stats::model.matrix(attr(mf,"terms"), data = mf)
  if(is.null(fe.spline)){
    fe.spline   <- list(df = 16, degree = 3, xi = 0.01)
  }
  if(is.null(se.spline)){
    se.spline   <- list(df = 8, degree = 3, slength = 100)
  }
  if(!is.null(fe)){
    if(is.null(Ts)){
      stop('Must specify grid for functional effect.')
    }
    cfe     <- construct_fe(fe, fe.spline, data, Ts)

    O           <- fe.O   <- cfe$O
    Th          <- fe.Th  <- cfe$Th
    Zf          <- fe.Zf  <- cfe$Zf
    data.names  <- fe.dn  <- names(stats::model.frame(formula = fe, data = data))
    Kfi         <- c(0, rep(fe.spline$df, length(O)))
    test.flag   <- rep('fe', length(fe.O))
  }
  if(!is.null(se)){
    cse     <- construct_se(se, se.spline, data)

    O                 <- se.O   <- cse$O
    Th                <- se.Th  <- cse$Th
    Zf                <- se.Zf  <- cse$Zf
    se.spline$sRange  <- cse$sRange
    data.names        <- se.dn  <- names(stats::model.frame(formula = se, data = data))
    Kfi               <- c(0, rep(se.spline$df, length(O)))
    test.flag         <- rep('se', length(se.O))
  }
  if(!is.null(fe) & !is.null(se)){
    O           <- c(se.O, fe.O)
    Th          <- c(se.Th, fe.Th)
    Zf          <- cbind(se.Zf, fe.Zf)
    data.names  <- c(se.dn, fe.dn)
    Kfi         <- c(0, rep(se.spline$df, length(se.O)), rep(fe.spline$df, length(fe.O)))
    test.flag   <- c(rep('se', length(se.O)), rep('fe', length(fe.O)))
  }
  n       <- length(Y) # number of subjects

  ### build design matrix ###
  XZ      <- cbind(Xc, Zf)

  ## prior values ##
  Ae    <- priors$Ae
  Be    <- priors$Be
  Af    <- priors$Af # edit defaults
  Bf    <- priors$Bf # edit defaults
  sigB  <- priors$sigB

  ## initialize algorithm ##
  Bqe         <- 1
  Bqf         <- rep(1, length = length(O))
  logp_delta	<- 0.5
  iter		    <- 0
  logp_prev	  <- 0

  ## set matrices ##
  XtX         <- t(Xc)%*%Xc
  CtC         <- t(XZ)%*%XZ
  p           <- ncol(Xc)
  Kf          <- ncol(Zf)
  MuMat       <- matrix(0, nrow = maxIter, ncol = p + Kf)
  SigmaMat    <- array(0, dim = c(p + Kf, p + Kf, maxIter))
  BqfMat      <- matrix(0, nrow = maxIter, ncol = length(O))
  lDelta      <- rep(logp_delta, maxIter)
  logpi       <- rep(0, maxIter)
  logpp       <- rep(0, maxIter)
  Sigq        <- vector('list', length = length(O) + 1)

  if(family$family == "gaussian"){
    CtY         <- t(XZ)%*%Y
    BqeVec      <- rep(0, maxIter)

    ## algorithm ##
    while(logp_delta > tol){
      iter	<- iter + 1

      ## update \Sigma_q ##
      Sigqb   <- Sigq[[1]]   <- (1/sigB)*diag(p)
      for(f in 1:length(O)){
        # Sigq[[f+1]] <- ((Af[f] + (1/2)*Kf/length(O))/Bqf[f])*O[[f]]
        Sigq[[f+1]] <- ((Af[f] + (1/2)*(Kfi[f + 1]))/Bqf[f])*O[[f]]
      }
      G       <- Matrix::bdiag(Sigq)
      Precq   <- as.matrix(((Ae + n/2)/Bqe)*CtC + G)
      Sigmaq  <- solve(Precq)

      ## update \mu_q ##
      Muq     <- (Ae + n/2)/(Bqe)*Sigmaq%*%CtY
      Ctf     <- (Ae + n/2)/(Bqe)*Sigmaq%*%t(XZ)

      ## update B_qe ##
      SSY     <- as.numeric(t(Y - XZ%*%Muq)%*%(Y - XZ%*%Muq))
      Bqe     <- Be + (1/2)*(SSY + sum(diag(CtC%*%Sigmaq)))

      ## update B_qf ##
      for(f in 1:length(O)){
        inds    <- (Kfi[f] + 1 + 1):(sum(Kfi[f:(f + 1)]) + 1)
        MltMf   <- as.numeric(t(Muq[inds])%*%Muq[inds])
        Bqf[f]  <- Bf[1] + (1/2)*(MltMf + sum(diag(Sigq[[f+1]])))
      }

      ## check criteria ##
      logp_i      <- eval_logp_gauss(Sigmaq = Sigmaq, sigB = sigB, Muq = Muq, p = p,
                                     Sigqb = Sigqb, Ae = Ae, Be = Be, n = n, Bqe = Bqe,
                                     Af = Af, Bf = Bf, Kfi = Kfi, Bqf = Bqf)

      logp_delta	<- abs(logp_i - logp_prev)

      ## track convergence ##
      lDelta[iter]  <- logp_delta
      logpi[iter]   <- logp_i
      logpp[iter]   <- logp_prev
      logp_prev	    <- logp_i

      ## store estimates ##
      MuMat[iter,]      <- as.numeric(Muq)
      SigmaMat[,,iter]  <- as.matrix(Sigmaq)
      BqeVec[iter]      <- Bqe
      BqfMat[iter,]     <- Bqf

      ## console update ##
      if(verbose){
        if(numbers::mod(iter, dots) == 0){
          cat('.')
        }

        if(numbers::mod(iter, up) == 0){
          cat(paste("\n",iter,"samples completed\n"))
        }
      }

      if(iter == maxIter){
        warning('Max number of iterations reached; increase max'); break
      }

    }
  } else{
    Muqa      <- rep(0, n)
    MuqaMat   <- matrix(0, nrow = maxIter, ncol = n)

    ## algorithm ##
    while(logp_delta > tol){
      iter	<- iter + 1

      ## update CtY
      CtY <- t(XZ)%*%Muqa

      ## update \Sigma_q ##
      Sigqb   <- Sigq[[1]]   <- (1/sigB)*diag(p)
      for(f in 1:length(O)){
        Sigq[[f+1]] <- ((Af[f] + (1/2)*(Kfi[f + 1]))/Bqf[f])*O[[f]]
      }
      G       <- Matrix::bdiag(Sigq)
      Precq   <- as.matrix(CtC + G)
      Sigmaq  <- solve(Precq)

      ## update \mu_q ##
      Muq     <- Sigmaq%*%CtY
      Ctf     <- Sigmaq%*%t(XZ)

      ## update B_qf ##
      for(f in 1:length(O)){
        inds    <- (Kfi[f] + 1 + 1):(sum(Kfi[f:(f + 1)]) + 1)
        MltMf   <- as.numeric(t(Muq[inds])%*%Muq[inds])
        Bqf[f]  <- Bf[1] + (1/2)*(MltMf + sum(diag(Sigq[[f+1]])))
      }

      ## update Muqa ##
      XZMuqa  <- XZ%*%Muq
      Muqa    <- XZMuqa + stats::dnorm(XZMuqa)/((stats::pnorm(XZMuqa)^Y)*((stats::pnorm(XZMuqa) - 1)^(1-Y)))

      ## check criteria ##
      logp_i      <- eval_logp_binpr(n = n, Y = Y, XZ = XZ, Muq = Muq, p = p,
                                     Sigqb = Sigmaq[1:p,1:p], Sigmaq = Sigmaq,
                                     Af = Af, Bf = Bf, Kfi = Kfi, Bqf = Bqf)

      logp_delta	<- abs(logp_i - logp_prev)

      ## track convergence ##
      lDelta[iter]  <- logp_delta
      logpi[iter]   <- logp_i
      logpp[iter]   <- logp_prev
      logp_prev	    <- logp_i

      ## store estimates ##
      MuMat[iter,]      <- as.numeric(Muq)
      MuqaMat[iter,]    <- as.numeric(Muqa)
      SigmaMat[,,iter]  <- as.matrix(Sigmaq)
      BqfMat[iter,]     <- Bqf

      ## console update ##
      if(verbose){
        if(numbers::mod(iter, dots) == 0){
          cat('.')
        }

        if(numbers::mod(iter, up) == 0){
          cat(paste("\n",iter,"samples completed\n"))
        }
      }

      if(iter == maxIter){
        warning('Max number of iterations reached; increase max'); break
      }

    }
  }

  ## algo outs ##
  algOut  <- list(maxIter = maxIter, tol = tol, totalIter = iter,
                  logpd = list(lDelta = lDelta, logp = logpi, prev_logp = logpp))

  ## outputs ##
  xcoef     <- matrix(MuMat[iter,1:p], p, 1)
  colnames(xcoef) <- "Estimate"
  rownames(xcoef) <- colnames(Xc)

  dmat      <- list(XZ = XZ, Xc = Xc, Zf = Zf, Ts = Ts)
  dn        <- vector('list', length = length(Th))
  for(j in 1:length(Th)){
    dn[[j]]   <- paste(c(names(mf)[1], data.names[j]), collapse = ' by ')
  }
  if(!is.null(fe)){
    fee       <- vector('list', length = length(fe.Th))
    for(j in 1:length(fe.Th)){
      inds      <- (Kfi[j] + 1 + 1):(sum(Kfi[j:(j + 1)]) + 1)
      fee[[j]]  <- as.numeric(Th[[j]]%*%MuMat[iter,inds])
    }
  } else{
    fee   <- NULL
  }
  if(!is.null(se)){
    see       <- vector('list', length = length(se.Th))
    for(j in 1:length(se.Th)){
      inds      <- (Kfi[j] + 1 + 1):(sum(Kfi[j:(j + 1)]) + 1)
      see[[j]]  <- as.numeric(Th[[j]]%*%MuMat[iter,inds])
    }
  } else{
    see   <- NULL
  }
  if(family$family == "gaussian"){
    BqeOut    <- BqeVec[iter]
    MuqaOut   <- NULL
  } else{
    BqeOut    <- NULL
    MuqaOut   <- MuqaMat[iter,]
  }
  out       <- list(VBEstimates = MuMat[iter,],
                    Muqa = MuqaOut,
                    fixedEffects = xcoef,
                    smoothEffect = see,
                    funcEffect = fee,
                    Ctf = Ctf, data.name = dn,
                    varComps = list(Sigmaq = SigmaMat[,,iter], Bqe = BqeOut,
                                    Bqf = BqfMat[iter,]),
                    model = list(fixed = fixed, fe = fe, se = se,
                                 data = data, p = p, Kf = Kf, Kfi = Kfi,
                                 family = family, test.flag = test.flag,
                                 dmat = dmat, Th = Th, O = O, Y = Y),
                    priors = priors, nfun = length(test.flag),
                    algOut = algOut,
                    se.spline = se.spline,
                    fe.spline = fe.spline)


  class(out) <- 'vam'
  return(out)

}
#' @export
plot.vam  <- function(x, which = 1:(x$nfun + 1),
                      ask = prod(graphics::par("mfcol")) < length(which) && grDevices::dev.interactive(),
                      add.points = TRUE, xlab = NULL, ylab = NULL,
                      ...){
  oldPar  <- graphics::par()
  oldAsk  <- oldPar$ask
  graphics::par(ask = ask)
  test.flag   <- x$model$test.flag
  if(!is.null(x$model$se)){
    mf      <- stats::model.frame(formula = x$model$fixed, data = x$model$data)
    mse     <- stats::as.formula(paste(deparse(x$model$se),"- 1"))
    Zl      <- stats::model.matrix(mse, data = x$model$data)
    se.dn   <- names(stats::model.frame(formula = x$model$se, data = x$model$data))
  }
  if(!is.null(x$model$fe)){
    fe.dn  <- names(stats::model.frame(formula = x$model$fe, data = x$model$data))
  }
  if(length(which) > 1){
    fe.counter  <- 0
    for(i in 1:x$nfun){
      if(test.flag[i] == 'se'){
        sRange    <- x$se.spline$sRange[[i]]
        if(is.null(xlab)){
          xlab  <- se.dn[i]
        }
        if(is.null(ylab)){
          ylab  <- names(mf)[1]
        }
        plot(sRange, x$model$family$linkinv(x$smoothEffect[[i]] + x$fixedEffects[1]), type = 'l',
             ylim = range(x$model$Y), ylab = ylab, xlab = xlab, ...)
        if(add.points){
          graphics::points(Zl[,i], x$model$Y)
        }
        xlab <- NULL
      } else{
        fe.counter <- fe.counter + 1
        if(is.null(xlab)){
          xlab  <- 'time'
        }
        if(is.null(ylab)){
          ylab  <- fe.dn[fe.counter]
        }
        plot(x$model$family$linkinv(x$funcEffect[[fe.counter]]), type = 'l', xlab = xlab, ylab = ylab, ...)
        ylab <- NULL
      }
    }
    plot(x$algOut$logpd$lDelta[1:x$algOut$totalIter], xlab = 'iteration',
         ylab = expression(paste(Delta, ' ', log(P))), type = 'b', pch = 16, ...)

    if(ask){
      graphics::par(ask = oldAsk)
    }
  } else{
    if(which > x$nfun){
      plot(x$algOut$logpd$lDelta[1:x$algOut$totalIter], xlab = 'iteration',
           ylab = expression(paste(Delta, ' ', log(P))), type = 'b', pch = 16, ...)
    } else{
      if(test.flag[which] == 'se'){
        sRange    <- x$se.spline$sRange[[which]]
        if(is.null(xlab)){
          xlab  <- se.dn[which]
        }
        if(is.null(ylab)){
          ylab  <- names(mf)[1]
        }
        plot(sRange, x$model$family$linkinv(x$smoothEffect[[which]] + x$fixedEffects[1]), type = 'l',
             ylim = range(x$model$Y), ylab = ylab, xlab = xlab, ...)
        if(add.points){
          graphics::points(Zl[,which], x$model$Y)
        }
      } else{
        if(is.null(xlab)){
          xlab  <- 'time'
        }
        if(is.null(ylab)){
          ylab  <- fe.dn[which]
        }
        plot(x$model$family$linkinv(x$funcEffect[[which]]), type = 'l', xlab = xlab, ylab = ylab, ...)
      }
    }
  }
}
#' @export
coef.vam  <- function(model){
  model$fixedEffects
}

