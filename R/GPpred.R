#' Conversion from a long dataframe to a sptm object
#'
#' This function allows you to calculate a Gaussian process kernel
#' @param xd Location of the GP prediction
#' @param x support points
#' @param y support values
#' @param param list of parameters, with the right naming convention
#' @param kernel list of kernel used to generate the GPs.
#' @keywords Gaussian process
#' @export
#' @examples
#' x = list(matrix(runif(5,0,10), ncol = 1))
#' y = runif(5,-2,2)
#' xd = list(matrix(seq(min(x[[1]])-1,max(x[[1]])+1,length.out = 100), ncol = 1))
#' param = list(list(q1=1, q2=2))
#' kernel = list(k_longterm)
#' t = GPpred(xd, x,y, param, kernel)
#'
#' xd = list(expand.grid(seq(0,10,length.out = 10),seq(0,10,length.out = 10)))
#' param = list(list(s1=2, s2=2))
#' x = list(matrix(runif(20,0,10),ncol = 2))
#' y = cos(rowSums(x[[1]]))
#' kernel = list(k_spatial_iso)
#' attr(kernel, "type") <- "spatial"
#' tt = GPpred(xd, x,y, param, kernel)
#' image(matrix(tt$mp, ncol = 10,byrow=TRUE))
GPpred = function(xd,x,y,param, kernel){

  y = matrix(y, ncol = 1)

    Kpx = matrix(0, nrow = nrow(xd[[1]]), ncol = nrow(x[[1]]))
    Kxx = matrix(0, nrow = nrow(x[[1]]), ncol = nrow(x[[1]]))
    Kpp = matrix(0, nrow = nrow(xd[[1]]), ncol = nrow(xd[[1]]))

    for(i in 1:length(kernel)){

      if(ncol(x[[i]])==1){
        Kpx = Kpx + outer(c(xd[[i]]),c(x[[i]]),kernel[[i]],param = param[[i]])
        Kxx = Kxx + outer(c(x[[i]]),c(x[[i]]),kernel[[i]],param = param[[i]])
        Kpp = Kpp + outer(c(xd[[i]]),c(xd[[i]]),kernel[[i]],param = param[[i]])
      }
      if(ncol(x[[i]])>1){
        xd1 = xd[[i]][rep(1:nrow(xd[[i]]), nrow(xd[[i]])),]
        xd2 = xd[[i]][rep(1:nrow(xd[[i]]), each = nrow(xd[[i]])),]
        x1 = x[[i]][rep(1:nrow(x[[i]]), nrow(x[[i]])),]
        x2 = x[[i]][rep(1:nrow(x[[i]]), each = nrow(x[[i]])),]

        xdd = cbind(xd1,xd2)
        xdx = cbind(xd[[i]][rep(1:nrow(xd[[i]]), nrow(x[[i]])),],x[[i]][rep(1:nrow(x[[i]]), each = nrow(xd[[i]])),])
        xx = cbind(x1,x2)
        Kpx = Kpx + matrix(apply(xdx,1, function(x) kernel[[i]](x[1:2],x[3:4],param = param[[i]])), ncol = nrow(x[[i]]))
        Kxx = Kxx +matrix(apply(xx,1, function(x) kernel[[i]](x[1:2],x[3:4],param = param[[i]])), ncol = nrow(x[[i]]))
        Kpp = Kpp + matrix(apply(xdd,1, function(x) kernel[[i]](x[1:2],x[3:4],param = param[[i]])), ncol = nrow(xd[[i]]))
      }

    }


  mu_p = Kpx %*% ginv(Kxx) %*% y
  sigma_p = Kpp - Kpx %*% ginv(Kxx) %*% t(Kpx)

  wd = list(mp = mu_p, sp = sigma_p)

  class(wd) <- "gppred"
  attr(wd, "type") <- attr(kernel, "type")

  return(wd)
}
