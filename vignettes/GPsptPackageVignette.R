## ------------------------------------------------------------------------
library(GPspt, quietly = TRUE)
x = seq(0,10,length.out = 100)
param = list(list(q1=2, q2=5))
nb = 10
kernel = list(k_longterm)
attr(kernel, "type") <- "temporal"
tt = GPgen(x, param, nb, kernel)
plot(tt$x, colMeans(tt$f), type = "l", col = "blue", ylim = range(tt$f), 
     xlab = "x", ylab = "f(x)")
for(i in 1:10){
  points(tt$x, tt$f[i,], type= "l", col = "lightgrey")
}

## ------------------------------------------------------------------------
set.seed(1)
x = list(matrix(runif(5,0,10), ncol = 1))
y = runif(5,-2,2)
xPred = list(matrix(seq(min(x[[1]])-1,max(x[[1]])+1,length.out = 100), ncol = 1))
param = list(list(q1=1, q2=2))
kernel = list(k_longterm)
tt = GPpred(xPred, x,y, param, kernel)
CI =  c(tt$mp -1.96*sqrt(diag(tt$sp)), rev(tt$mp +1.96*sqrt(diag(tt$sp))))
plot(xPred[[1]], tt$mp, type = "l", col = "blue", xlab = "x", ylab = "f(x)", ylim = range(CI))
points(x[[1]],y, col="red", pch = 3)
polygon(c(xPred[[1]], rev(xPred[[1]])),CI, col = "lightgrey", density = 30)

## ------------------------------------------------------------------------
kernelList = list(k_longterm)
attr(kernelList, "name") <- c("long term")
attr(kernelList, "parameters") <- list(c("q1","q2"))
attr(kernelList, "type") <- c("temporal")
thetaInit = c(log(2),log(1))
x  = matrix(seq(0,10,1), ncol = 1)
y = cos(x)

xList = list(x)
res = optifix(thetaInit, fixed = c(FALSE,FALSE),fn = wrap_ll, y=y, kernel = kernelList, 
              x=xList,  control = list(fnscale = -1, maxit = 10000), method = "Nelder-Mead")
param = parToList(res$fullpar, kernelList)

xPred = list(matrix(seq(min(x)-1,max(x)+1,length.out = 100), ncol = 1))
tt = GPpred(xPred, xList,y, param, kernelList)
CI =  c(tt$mp -1.96*sqrt(diag(tt$sp)), rev(tt$mp +1.96*sqrt(diag(tt$sp))))
plot(xPred[[1]], tt$mp, type = "l", col = "blue", xlab = "x", ylab = "f(x)", ylim = range(CI))
points(xList[[1]],y, col="red", pch = 3)
polygon(c(xPred[[1]], rev(xPred[[1]])),CI, col = "lightgrey", density = 30)

## ------------------------------------------------------------------------
kernelList = list(k_spatial_iso)
attr(kernelList, "name") <- c("spatial")
attr(kernelList, "parameters") <- list(c("s1", "s2"))
attr(kernelList, "type") <- c("spatial")
thetaInit = c(log(1), log(5), log(0.1))
s = expand.grid(seq(1,10,length.out=5), seq(1,10,length.out = 5))
y = exp(-1/25*((s[,1]-5)^2 + (s[,2]-5)^2))
xList = list(s)
res = optifix(thetaInit, fixed = c(FALSE, FALSE,FALSE),fn = wrap_ll, y=y, kernel = kernelList, x=xList,  control = list(fnscale = -1, maxit = 10000), method = "Nelder-Mead")
param = parToList(res$fullpar, kernelList)
#param = parToList(thetaInit, kernelList)

## ------------------------------------------------------------------------
sPred = expand.grid(seq(1,10,length.out=20), seq(1,10,length.out = 20))

xPred = list(sPred)
tt = GPpred(xPred, xList,y, param, kernelList)

xPred2 =  list(s)
tt2 = GPpred(xPred2, xList,y, param, kernelList)

p <- persp(x =seq(1,10,length.out=20) , y=seq(1,10,length.out=20) , z=matrix(tt$mp, ncol=20), phi = 30, theta=40, xlab = "x", ylab = "y", zlab = "z",col="lightblue",expand = 0.5,shade = 0.2)
obs <- trans3d(xList[[1]][,1], xList[[1]][,2],y,p)
pred <- trans3d(xPred2[[1]][,1], xPred2[[1]][,2],tt2$mp,p)
points(obs, col="red",pch=16)
points(pred, col="blue",pch=16, cex=0.5)
segments(obs$x, obs$y, pred$x, pred$y)

