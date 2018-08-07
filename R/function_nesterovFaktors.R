#' Nesterov Faktor 2/(k+1)
#'
#' The FISTA algorithm extrapolates/correctes after the gradient step. The "nesterov_faktor" function
#' returns how far the FISTA will exptrapolate in a given step "k"
#'
#' @param k integer, iteration step of the FISTA algorithm
#'
#' @return
#'
#' @export
#'
#' @examples
#' nesterov_faktor(2)
nesterov_faktor <- function(k){
  return(2/(k+1))
}

