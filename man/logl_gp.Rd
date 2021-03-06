% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{logl_gp}
\alias{logl_gp}
\title{GP log-likelihood function}
\usage{
logl_gp(paramList, y, kernelList, xList)
}
\arguments{
\item{paramList}{list of parameters}

\item{y}{observations}

\item{kernelList}{list of kernels}

\item{xList}{list of features (time, location, features)}
}
\description{
GP log-likelihood function
}
\examples{
kernelList = list(k_seasonal, k_spatial_iso)
paramList =  list(list(q3=1,q4=1,qs=1, f=1), list(s1 = 1, s2 = 5), sigma_noise = 0.2)
t  = matrix(rep(seq(0,10,1),11), ncol = 1)
s = expand.grid(seq(0,10,1), seq(0,10,1))
y = cos(rowSums(s) * t)
xList = list(t,s)
ll = logl_gp(paramList, y, kernelList,xList)
}
\keyword{Gaussian}
\keyword{likelihood}
\keyword{process}
