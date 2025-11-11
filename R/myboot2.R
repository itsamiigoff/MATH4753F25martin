#' myboot2
#'
#' @param iter number of bootstrap replications
#' @param x sample data
#' @param fun statistic to compute
#' @param alpha significance level (1 - alpha is ci)
#' @param cx the plot's text size
#' @param ... optional plotting parameters
#' @importFrom graphics segments
#' @importFrom stats quantile
#'
#' @returns a list with ci, the function used, and the original data
#' @export
#'
#' @examples
#' myboot2(x=ddt$DDT)
myboot2<-function(iter=10000,x,fun="mean",alpha=0.05,cx=1.5,...){
  # sample size
  n=length(x)

  y=sample(x,n*iter,replace=TRUE) # Line A
  rs.mat=matrix(y,nrow=n,ncol=iter,byrow=TRUE)

  # xstat is a vector and will have iter values in it
  xstat=apply(rs.mat,2,fun)

  # Nice way to form a confidence
  ci=quantile(xstat,c(alpha/2,1-alpha/2))# Line B interval

  # A histogram follows
  # The object para will contain the parameters used to make the histogram
  para=hist(xstat,freq=FALSE,las=1,
            main=paste("Histogram of Bootstrap sample statistics", "\n", "alpha=", alpha, "iter=", iter, sep=""),
            ...)

  #mat will be a matrix that contains the data, this is done so that I can use apply()
  mat=matrix(x,nrow=length(x),ncol=1,byrow=TRUE)

  #pte is the point estimate
  #This uses whatever fun is
  pte=apply(mat,2,fun)

  # Vertical line
  abline(v=pte,lwd=3,col="Black")

  # Make the segment for the ci
  segments(ci[1],0,ci[2],0,lwd=4)
  text(ci[1],0,paste("(",round(ci[1],2),sep=""),col="Red",cex=cx)
  text(ci[2],0,paste(round(ci[2],2),")",sep=""),col="Red",cex=cx)

  # plot the point estimate 1/2 way up the density
  text(pte,max(para$density)/2,round(pte,2),cex=cx)

  # Some output to use if necessary
  invisible(list(ci=ci,fun=fun,x=x))
}
