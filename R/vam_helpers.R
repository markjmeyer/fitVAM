construct_P  <- function(D, xi, Theta){
  diff0   <- diag(1, D, D)
  diff2   <- matrix(rep(c(1, -2, 1, rep(0, D - 2)), D - 2)[1:((D - 2) * D)], D - 2, D,
                    byrow = TRUE)
  P0      <- t(Theta) %*% t(diff0) %*% diff0 %*% Theta
  P2      <- t(Theta) %*% t(diff2) %*% diff2 %*% Theta
  P.mat   <- xi * P0 + (1 - xi) * P2

  return(P.mat)
}

construct_fe  <- function(fe, ps.spline, data, Ts){
  mfe     <- stats::as.formula(paste(deparse(fe),"- 1"))
  Zl      <- stats::model.matrix(mfe, data = data)
  dfs     <- ps.spline$df
  degs    <- ps.spline$degree
  xi      <- ps.spline$xi
  Ks      <- vector(length = length(Ts))
  Kf      <- length(Ts)*dfs
  Th      <- O <-   vector('list', length = length(Ts))
  Zf      <- matrix(NA, nrow(Zl), Kf)
  Tsi     <- c(0, unlist(Ts))
  Kfi     <- c(0, rep(dfs, length(Ts)))
  for(i in 1:length(Ts)){
    Ks[i]     <- Ts[[i]]
    x         <- 1:Ks[i]
    Th[[i]]   <- splines::bs(x, df = dfs, intercept = FALSE, degree = degs)
    O[[i]]    <- construct_P(Ks[i], xi, Th[[i]])

    Zf[,(Kfi[i] + 1):(sum(Kfi[i:(i+1)]))]      <- Zl[,(Tsi[i] + 1):(sum(Tsi[i:(i+1)]))]%*%Th[[i]]
  }

  out     <- list(Zf = Zf, Th = Th, O = O)

  return(out)
}

construct_se  <- function(se, se.spline, data){
  mse     <- stats::as.formula(paste(deparse(se),"- 1"))
  Zl      <- stats::model.matrix(mse, data = data)
  dfs     <- se.spline$df
  Kf      <- ncol(Zl)*dfs
  degs    <- se.spline$degree
  slength  <- se.spline$slength
  Zf      <- matrix(NA, nrow(Zl), Kf)
  Kfi     <- c(0, rep(dfs, ncol(Zl)))
  Th      <- O <- sRange  <-  vector('list', length = ncol(Zl))
  for(i in 1:ncol(Zl)){
    Zf[,(Kfi[i] + 1):(sum(Kfi[i:(i+1)]))] <- Zfi <- splines::bs(as.numeric(Zl[,i]), df = dfs, intercept = FALSE, degree = degs)
    Ks      <- ncol(Zfi)
    Db      <- diff(diag(Ks))

    O[[i]]        <- t(Db)%*%Db
    sRange[[i]]   <- seq(range(Zl[,i])[1], range(Zl[,i])[2], length = slength)
    Th[[i]]       <- splines::bs(sRange[[i]], df = dfs, intercept = FALSE, degree = degs)
  }

  out     <- list(Zf = Zf, Th = Th, O = O, sRange = sRange)

  return(out)
}

eval_logp_gauss  <- function(Sigmaq, sigB, Muq, p, Sigqb, Ae, Be, n, Bqe, Af, Bf, Kfi, Bqf){
  if(n/2 < 150){
    out   <- (1/2)*(p + sum(Kfi[-1])) - (n/2)*log(2*pi) - (p/2)*log(sigB) +
      (1/2)*log(det(Sigmaq)) - (1/(2*sigB))*(t(Muq[1:p])%*%Muq[1:p] + sum(diag(Sigqb))) -
      Ae*log(Be) - (Ae + n/2)*log(Bqe) + log(gamma(Ae + n/2)) - log(gamma(Ae)) +
      sum(Af*log(Bf) - (Af + Kfi[-1]/2)*log(Bqf) + log(gamma(Af + Kfi[-1]/2)) - log(gamma(Af)))
  } else{
    out   <- (1/2)*(p + sum(Kfi[-1])) - (n/2)*log(2*pi) - (p/2)*log(sigB) +
      (1/2)*log(det(Sigmaq)) - (1/(2*sigB))*(t(Muq[1:p])%*%Muq[1:p] + sum(diag(Sigqb))) -
      Ae*log(Be) - (Ae + n/2)*log(Bqe) - log(gamma(Ae)) +
      sum(Af*log(Bf) - (Af + Kfi[-1]/2)*log(Bqf) + log(gamma(Af + Kfi[-1]/2)) - log(gamma(Af)))
  }

  return(as.numeric(out))
}

eval_logp_binpr <- function(n, Y, XZ, Muq, p, Sigqb, Sigmaq, Af, Bf, Kfi, Bqf){
  ones  <- rep(1, n)
  out   <- t(Y)%*%log(stats::pnorm(XZ%*%Muq) + 1) +
    t(ones - Y)%*%log(ones - stats::pnorm(XZ%*%Muq) + 1) -
    (1/(2))*(t(Muq[1:p])%*%Muq[1:p] + sum(diag(Sigqb))) -
    (1/2)*log(det(Sigmaq)) +
    sum(Af*log(Bf) - (Af + Kfi[-1]/2)*log(Bqf) + log(gamma(Af + Kfi[-1]/2)) - log(gamma(Af)))
}
