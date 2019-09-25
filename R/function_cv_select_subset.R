#' selects a specific sample from a list of entries
#'
#' @param list.entry list of entries (matrix or list)
#' @param samples samples, serving as an index to list.entry.
#'
#' @return specific entry of list.entry
#'
select.fun <- function(list.entry, samples) {
  # internal training samples selection function:
  if (is.matrix(list.entry)) {
    return(list.entry[, samples])
  }
  if (!is.null(names(list.entry))) {
    return(list.entry[samples])
  }
  return(list.entry)
}
