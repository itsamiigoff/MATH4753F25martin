#' ddt data
#'
#' A dataset containing DDT concentrations in fish collected from several rivers.
#' Each record includes river location, species, length, weight, and measured DDT level.
#'
#' @format A data frame with 120 rows and 6 variables:
#' \describe{
#'   \item{RIVER}{River name or sampling site (e.g., FCM, LCM, SCM, TRM)}
#'   \item{MILE}{River mile (numeric) indicating sampling location}
#'   \item{SPECIES}{Fish species (e.g., CCATFISH, SMBUFFALO, LMBASS)}
#'   \item{LENGTH}{Fish length in centimeters}
#'   \item{WEIGHT}{Fish weight in grams}
#'   \item{DDT}{Concentration of DDT (µg/g) in fish tissue}
#' }
#'
#' @source Example dataset commonly used in environmental statistics labs (adapted from U.S. Geological Survey sample data)
"ddt"
