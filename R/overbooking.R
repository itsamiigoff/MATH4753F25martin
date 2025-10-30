#' Overbooking
#'
#' @param N number of seats on the flight
#' @param gamma probability of overbooking - the "pain" factor
#' @param p probability a passenger shows up
#' @importFrom graphics abline lines title
#' @importFrom stats pbinom qbinom qnorm uniroot
#'
#' @returns list containing optimal number of seats sold via discrete and
#' continuous models, plots the two models as well
#' @export
#'
#' @examples
#' ntickets(N=400,gamma = 0.02, p = 0.95)
ntickets = function(N, gamma, p) {
  # Optimal n for discrete using qbinom
  nd = 2 * N - qbinom(1 - gamma, size = N, prob = p)

  # Optimal n for continuous using qnorm
  nc = 2 * N - qnorm(1 - gamma, mean = N * p, sd = sqrt(N * p * (1 - p)))

  # Discrete objective function f(n)
  discrete = function(n) {
    1 - gamma - pbinom(N, round(n), p)
  }

  # Continuous objective function f(n)
  continuous = function(n) {
    1 - gamma - pnorm(N + 0.5, n * p, sqrt(n * p * (1 - p)))
  }

  # Evaluate discrete objective
  n_vals = seq(N, N + 20, by = 1) # Possible number of tickets to test
  obj_vals = sapply(n_vals, discrete) # Apply discrete to each value in n_vals

  # Plot points and connect with a line
  plot(n_vals, obj_vals, type = "p", pch = 19,
       xlab = "n", ylab = "Objective", ylim = c(0, 1))
  lines(n_vals, obj_vals, col = "black")

  # Find approximate optimal n for discrete
  unirootDis = uniroot(f = discrete, interval = c(N - 5, N + 25))$root

  # Add reference lines
  abline(v = round(unirootDis), h = 0, col = "red", lwd = 2)

  # Add title
  title(main = paste0(
    "Objective Vs n to find optimal tickets sold\n(",
    round(unirootDis), ") gamma= ", gamma, " N=", N, " discrete"
  ))


  # Plot continuous
  dplot = curve(expr = continuous, from = N, to = N + 20, col = "black",
                xlab="n",ylab="Objective")

  # Find optimal n for continuous
  unirootCont = uniroot(f = continuous, interval = c(N - 5, N + 25))$root

  # Add reference lines
  abline(v = unirootCont, h = 0, lty = "solid", col = "blue")

  # Add title
  title(main=paste0("Objective Vs n to find optimal tickets sold \n(",
                    nc, "), gamma = ",gamma, ", N = ", N, ", continuous"))

  # Return list of results
  result = list(nd = nd, nc = nc, N = N, p = p, gamma = gamma)
  return(result)
}
