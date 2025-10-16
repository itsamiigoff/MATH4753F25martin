#' myncurve
#'
#' @param mu the mean of the function
#' @param sigma the standard deviation of the function
#' @param a the x value to compute the probability
#' @importFrom graphics curve polygon text
#' @importFrom stats dnorm pnorm
#'
#' @returns the probability P(X≤a) and the graph
#' @export
#'
#' @examples
#' myncurve(mu=10,sigma=5, a=6)
myncurve = function(mu, sigma, a){
  curve(dnorm(x,mean=mu,sd=sigma), xlim = c(mu-3*sigma, mu + 3*sigma))

  xcurve = seq(mu - 4 * sigma, a, length = 1000)
  ycurve = dnorm(xcurve, mean = mu, sd = sigma)

  polygon(c(mu-4*sigma,xcurve,a), c(0,ycurve,0), col = "#00FF00")
  area = pnorm(a, mean = mu, sd = sigma)
  area = round(area, 4)

  text(x = -2.5, y = 0.05, paste0("Area = ", area))
  return(list(mean = mu, sd = sigma, a = a, area = area))
}
