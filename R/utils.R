#' GP log-likelihood function
#'
#' @param paramList list of parameters
#' @param y observations
#' @param kernelList list of kernels
#' @param xList list of features (time, location, features)
#' @keywords Gaussian process likelihood
#' @export
#' @examples
#' kernelList = list(k_seasonal, k_spatial_iso)
#' paramList =  list(list(q3=1,q4=1,qs=1, f=1), list(s1 = 1, s2 = 5), sigma_noise = 0.2)
#' t  = matrix(rep(seq(0,10,1),11), ncol = 1)
#' s = expand.grid(seq(0,10,1), seq(0,10,1))
#' y = cos(rowSums(s) * t)
#' xList = list(t,s)
#' ll = logl_gp(paramList, y, kernelList,xList)
logl_gp = function(paramList, y, kernelList,xList){
  y = c(y)

  Ky = K_wrap(paramList, y, kernelList,xList)

  ll = -1/2 * t(y) %*% ginv(Ky) %*% y - 1/2 * as.numeric(determinant(Ky, logarithm = TRUE)$modulus) - length(y) /2 * log(2*pi)
  return(ll)
}

#' Wrapper for kernel calculation
#'
#' @param paramList list of parameters
#' @param y observations
#' @param kernelList list of kernels
#' @param xList list of features (time, location, features)
#' @keywords Gaussian process likelihood
#' @export
#' @examples
#' kernelList = list(k_spatial_iso)
#' paramList =  list(list(s1 = 1, s2 = 2), sigma_noise = 0.2)
#' s = expand.grid(seq(0,10,1), seq(0,10,1))
#' y = exp( -1/5 * ((s[,1]-5)^2+(s[,2]-5)^2))
#' xList = list(s)
#' K = K_wrap(paramList, y, kernelList,xList)
#' persp(matrix(K[,3], ncol=11))
K_wrap = function(paramList, y, kernelList,xList){
  y = c(y)
  Kf = matrix(0, ncol = length(y), nrow = length(y))
  for(i in 1:length(kernelList)){
    x = xList[[i]]
    if(ncol(x) == 1){
      Ktemp = outer(c(x),c(x), kernelList[[i]], param = paramList[[i]])
    }else{
      xd1 = x[rep(1:nrow(x), nrow(x)),]
      xd2 = x[rep(1:nrow(x), each = nrow(x)),]
      xdd = cbind(xd1,xd2)
      Ktemp = matrix(apply(xdd,1, function(x) kernelList[[i]](x[1:2],x[3:4],param = paramList[[i]])), ncol = nrow(x))
    }
    Kf = Kf + Ktemp
  }
  Ky =  Kf + paramList$sigma_noise * diag(length(y))

  return(Ky)
}

#' Wrapper for GP likelihood maximisation
#'
#' @param theta parameters to be estimated
#' @param y observations
#' @param kernelList list of kernels
#' @param xList list of features
#' @keywords Gaussian process likelihood wrapper
#' @export
#' @examples
#' kernelList = list(k_longterm, k_seasonal, k_spatial)
#' attr(kernelList, "name") <- c("long term", "seasonal", "spatial")
#' attr(kernelList, "parameters") <- list(c("q1","q2"), c("q3","q4","qs", "f"), c("s1", "s2"))
#' theta = c(log(2),log(1),log(1),log(1),log(1),log(1),log(1), log(1), log(0.2))
#' t  = matrix(rep(seq(0,10,1),11), ncol = 1)
#' s = expand.grid(seq(0,10,length.out=11), seq(0,10,length.out = 11))
#' y = cos(rowSums(s) * t)
#' xList = list(t,t,s)
#' ll = wrap_ll(theta, y, kernelList,xList)
wrap_ll = function(theta,y,kernelList,xList){
  idx.p = 0
  param = NULL
  for(i in 1:length(kernelList)){
    paramTemp=NULL
    par = attr(kernelList, "parameters")[[i]]
    idx.p = max(idx.p) + seq_len(length(par))
    for(j in 1:length(idx.p)){
      assign(par[j], exp(theta[idx.p][j]))
    }
    eval(parse(text=paste0("paramTemp = c(paramTemp,c(list(",par,"=",par,")))")))
    param = c(param, list(paramTemp))
  }

  paramList = c(param, list(sigma_noise = exp(theta[length(theta)])))
  ll = logl_gp(paramList, y, kernelList = kernelList, xList=xList)
  return(ll)
}

#' Optimise with fixed parameters
#' @param par paramters
#' @param fixed boolean
#' @param fn function
#' @param gr gradient of the function
#' @param method optimization method
#' @param lower lower bounds.
#' @param upper upper bounds.
#' @param control list of controls.
#' @param hessian boolean. Return the Hessian.
#' @param ... additional parameters related to the fn and gr functions.
#' @keywords optimization with fixed values.
#' @export
#' @examples
#' fn = function(t){(2 - t[1] + t[2])^2 + 2*(3 - t[3])^4}
#' th = c(1,1,1)
#' opt = optifix(th,fixed = c(TRUE, FALSE, FALSE), fn = fn, method = "Nelder-Mead")
optifix <- function(par, fixed, fn, gr = NULL, ..., method = c("Nelder-Mead", "BFGS", "CG", "L-BFGS-B", "SANN"), lower = -Inf, upper = Inf, control = list(), hessian = FALSE){
  force(fn)
  force(fixed)
  .npar=length(par)
  .fixValues = par[fixed]
  .parStart = par[!fixed]
  .fn <- function(par,...){
    .par = rep(NA,sum(!fixed))
    .par[!fixed] = par
    .par[fixed] = .fixValues
    fn(.par,...) }
  if(!is.null(gr)){
    .gr <- function(par,...){
      .gpar = rep(NA,sum(!fixed))
      .gpar[!fixed] = par
      .gpar[fixed] = .fixValues
      gr(.gpar,...)[!fixed] }
  }else{ .gr <- NULL }

  .opt = optim(.parStart,.fn,.gr,...,method=method,lower=lower,control=control,hessian=hessian)
  .opt$fullpars = rep(NA,sum(!fixed))
  .opt$fullpars[fixed]=.fixValues
  .opt$fullpars[!fixed]=.opt$par
  .opt$fixed = fixed
  return(.opt)
}

#' Converts theta to hyperparameters
#'
#' @param theta parameters to be estimated
#' @param kernelList list of kernels
#' @keywords Gaussian process wrapper
#' @export
#' @examples
#' kernelList = list(k_longterm, k_seasonal, k_spatial)
#' attr(kernelList, "name") <- c("longterm", "seasonal", "spatial")
#' attr(kernelList, "parameters") <- list(c("q1","q2"), c("q3","q4","qs", "f"), c("s1", "s2"))
#' theta = c(log(1),log(1),log(1),log(1),log(1),log(1),log(1), log(1), log(0.2))
#' res = parToList(theta, kernelList)
parToList = function(theta, kernelList){
  idx.p = 0
  param = NULL
  for(i in 1:length(kernelList)){
    paramTemp=NULL
    par = attr(kernelList, "parameters")[[i]]
    idx.p = max(idx.p) + seq_len(length(par))
    for(j in 1:length(idx.p)){
      assign(par[j], exp(theta[idx.p][j]))
    }
    eval(parse(text=paste0("paramTemp = c(paramTemp,c(list(",par,"=",par,")))")))
    param = c(param, list(paramTemp))
  }

  param = c(param, list(sigma_noise = exp(theta[length(theta)])))
  names(param)[1:length(kernelList)] <- attr(kernelList, "name")
  return(param)
}
