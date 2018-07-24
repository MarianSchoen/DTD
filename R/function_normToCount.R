#' normalizeToCount
#' 
#' 'normalizeToCount' takes a numeric matrix mat and normalizes each column (=> sample) to fix number of counts
#'
#' @param mat numeric matrix, with samples as columns and features as columns
#' @param count float 
#'
#' @return : "ret", matrix with same dimension as "mat"
#' @export
#'
#' @examples 
#' someMatrix <- matrix(rnorm(1000*5), ncol=5, nrow=1000)
#' # each sample (=column) has different sum: 
#' apply(someMatrix, 2, sum)
#' 
#' normalizedMatrix <- normalizeToCount(someMatrix)
#' 
#' # check: 
#' apply(normalizedMatrix, 2, sum)
normalizeToCount <- function(mat, count=1e6){
  # get all sum over all columns: 
  sums <- colSums(mat)
  
  # for each column ...
  for(l1 in 1:ncol(mat)){
    # ... norm each column by its sum, then multiply each entry with "count"
    mat[,l1] <- count * mat[,l1]/sums[l1]
  }
  return(mat)
}
