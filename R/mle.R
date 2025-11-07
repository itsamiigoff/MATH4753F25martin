#' SHINY MLE
#'
#' @returns shiny app
#' @export
#'
#' @examples
#' \dontrun{shinymle()}
shinymle = function(){
  shiny::runApp(system.file("SHINY", "MLEapp", package = "MATH4753F25martin"),
  launch.browser = TRUE)
}
