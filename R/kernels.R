#' Long term tempporal kernel
#'
#' @param t Values for which the kernel is calculated
#' @param tp Values for which the kernel is calculated
#' @param param list of parameters specific to the kernel
#' @keywords Gaussian process kernel
#' @export
#' @examples
#' t = seq(0,1,length.out = 10)
#' param = list(q1=1, q2=1)
#' k = outer(t,t, k_longterm, param = param)
k_longterm = function(t,tp, param = list(q1=1,q2=1)){
  param$q1^2 * exp(-(t-tp)^2/ param$q2^2)
}
#' Seasonal tempporal kernel
#'
#' @param t Values for which the kernel is calculated
#' @param tp Values for which the kernel is calculated
#' @param param list of parameters specific to the kernel
#' @keywords Gaussian process kernel
#' @export
#' @examples
#' t = seq(0,1,length.out = 10)
#' param = list(q3=1, q4=1, qs = 1, f = 1)
#' k = outer(t,t, k_seasonal, param = param)
k_seasonal = function(t,tp,param = list(q3=1,q4=1,qs=1,f=1)){
  param$q3^2 * exp( -2*sin(3.14*param$f*(t-tp))^2 / param$qs^2) * exp(-1/2*(t-tp)^2 / param$q4^2)
}
#' Short term tempporal kernel
#'
#' @param t Values for which the kernel is calculated
#' @param tp Values for which the kernel is calculated
#' @param param list of parameters specific to the kernel
#' @keywords Gaussian process kernel
#' @export
#' @examples
#' t = seq(0,1,length.out = 10)
#' param = list(t6=1, t7=1 ,t8 = 1)
#' k = outer(t,t, k_shortterm, param = param)
k_shortterm = function(t,tp,param = list(t6=1,t7=1,t8=1)){
  param$t6^2* (1 + (t-tp)^2 / (2 * param$t7 * param$t8)) - param$t8
}
#' Spatial isotropic kernel
#'
#' @param s Locations for which the kernel is calculated
#' @param sp Locations for which the kernel is calculated
#' @param param list of parameters specific to the kernel
#' @keywords Gaussian process kernel
#' @export
#' @examples
#' x = seq(0,10,length.out = 4)
#' xd = expand.grid(x,x)
#' xd1 = xd[rep(1:nrow(xd), nrow(xd)),]
#' xd2 = xd[rep(1:nrow(xd), each = nrow(xd)),]
#' xdd = cbind(xd1,xd2)
#' param = list(s1 = 1, s2 = 10)
#' k = matrix(apply(xdd,1, function(x) k_spatial_iso(x[1:2],x[3:4],param = param)), ncol = nrow(xd))
k_spatial_iso= function(s,sp,param = list(s1=1,s2=1)){
  param$s1^2 * exp(-(t(s-sp) %*% (s-sp))/ param$s2^2)
}
#' Spatial Mattern kernel
#'
#' @param s Locations for which the kernel is calculated
#' @param sp Locations for which the kernel is calculated
#' @param param list of parameters specific to the kernel
#' @keywords Gaussian process kernel
#' @export
#' @examples
#' x = seq(0,10,length.out = 4)
#' xd = expand.grid(x,x)
#' xd1 = xd[rep(1:nrow(xd), nrow(xd)),]
#' xd2 = xd[rep(1:nrow(xd), each = nrow(xd)),]
#' xdd = cbind(xd1,xd2)
#' param = list(m1 = 1, m2 = 10)
#' k = matrix(apply(xdd,1, function(x) k_spatial_matern(x[1:2],x[3:4],param = param)), ncol = nrow(xd))
k_spatial_matern= function(s,sp,param = list(m1=1,m2=1)){
  d = sqrt(t(s-sp) %*% (s-sp))
  param$m1^2 * (1 + sqrt(3) *d / param$m2) * exp(- sqrt(3) *d / param$m2)
}
#' Features isotropic kernel
#'
#' @param f Features for which the kernel is calculated
#' @param fp Features for which the kernel is calculated
#' @param param list of parameters specific to the kernel
#' @keywords Gaussian process kernel
#' @export
#' @examples
#' f = matrix(rnorm(100,0,5), ncol = 5)
#' fd1 = f[rep(1:nrow(f), nrow(f)),]
#' fd2 = f[rep(1:nrow(f), each = nrow(f)),]
#' fdd = cbind(fd1,fd2)
#' param = list(s1 = 1, s2 = 10)
#' kern = apply(fdd,1, function(x) k_feature(x[1:ncol(f)],x[(ncol(f)+1):ncol(fdd)],param = param))
#' k = matrix(kern, ncol = nrow(f))
k_feature = function(f,fp, param = list(f1=1,f2=1)){
  param$f1^2 * exp(-(t(f-fp) %*% (f-fp))/ param$f2^2)
}
