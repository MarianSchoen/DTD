#' returns mean of the individual test results, setting unconverged or unavailable results to Inf.
#'
#' @param lambda.list output of a DTD cross validation run (list of lists)
#'
#' @return list of means.
#'
pick.mean.test.results.function <- function(lambda.list){
  tmp <- lapply(lambda.list, function(each.fold){
    if("cor.test" %in% names(each.fold)){
      return(each.fold$cor.test)
    }else{
      return(NA)
    }
  })
  test.vec <- mean(unlist(tmp, use.names = FALSE), na.rm = TRUE)
  return(test.vec)
}
