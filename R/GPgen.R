#' Conversion from a long dataframe to a sptm object
#'
#' This function allows you to calculate a Gaussian process kernel
#' @param xd Location of the GP prediction
#' @param param list of parameters, with the right naming convention
#' @param nb Number of GP to be generated
#' @param kernel list of kernel used to generate the GPs.
#' @keywords Gaussian process
#' @export
#' @examples
#' xd = expand.grid(seq(0,10,length.out = 10),seq(0,10,length.out = 10))
#' param = list(list(s1=2, s2=5))
#' nb = 1
#' kernel = list(k_spatial)
#' attr(kernel, "type") <- "spatial"
#' tt = GPgen(xd, param, nb, kernel)
#' image(matrix(tt$f, ncol = 10,byrow=TRUE))
GPgen = function(xd, param, nb, kernel){
  if(attr(kernel, "type") == "temporal"){
    sigma_pd = matrix(0, ncol = length(xd), nrow = length(xd))
    for(i in 1:length(kernel)){
      sigma_pd = sigma_pd + outer(xd,xd, kernel[[i]], param = param[[i]])
    }
    wd = list(x = xd, f = mvrnorm(nb,mu=rep(0, nrow(sigma_pd)),sigma_pd))
  }
  if(attr(kernel, "type") == "spatial"){
    sigma_pd = matrix(0, ncol = nrow(xd), nrow = nrow(xd))
    xd1 = xd[rep(1:nrow(xd), nrow(xd)),]
    xd2 = xd[rep(1:nrow(xd), each = nrow(xd)),]
    xdd = cbind(xd1,xd2)
    for(i in 1:length(kernel)){
      sigma_pd = sigma_pd + matrix(apply(xdd,1, function(x) kernel[[i]](x[1:2],x[3:4],param = param[[i]])), ncol = nrow(xd))
    }
    wd = list(x = xd, f = mvrnorm(nb,mu=rep(0, nrow(sigma_pd)),sigma_pd))
  }
  class(wd) <- "gpgen"
  attr(wd, "type") <- attr(kernel, "type")
  return(wd)
}
