#' Birthday
#'
#' @param x number of people in the room
#'
#' @returns The probability that two or more people in a class of size x have the same birthday.
#' @export
#'
#' @examples
#' birthday(20:25)
birthday <- function(x){
  1 - exp(lchoose(365,x) + lfactorial(x) - x*log(365))
}
