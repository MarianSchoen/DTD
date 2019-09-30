#' test_integer
#' (All test functions are used for safety checks in the DTD package)
#'
#' @param test.value value to be tested for integer
#' @param output.info [1] function that calls, [2] name of value to test
#' @param min min is tested '<'
#' @param max max is tested '>'
#'
#' @return TRUE, if no error is detected, stops with error otherwise
test_integer <- function(test.value,
                         output.info,
                         min,
                         max){

  error.message <- paste0("In ", output.info[1], ": ", "'", output.info[2], "'")

  if(!is.numeric(test.value) || length(test.value) != 1){
    error.message <- paste0(error.message, " is not a single integer")
    stop(error.message, call. = FALSE)
  }
  if(round(test.value) != test.value ){
    error.message <- paste0(error.message, " is not an integer")
    stop(error.message, call. = FALSE)
  }
  if(test.value < min){
    error.message <- paste0(error.message, " is below minimal value")
    stop(error.message, call. = FALSE)
  }
  if(test.value > max){
    error.message <- paste0(error.message, " is above maximal value")
    stop(error.message, call. = FALSE)
  }
  return(TRUE)
}


#' test_tweak_vec
#' (All test functions are used for safety checks in the DTD package)
#' @param tweak.vec vector which will be tested for numeric and NAs
#' @param output.info [1] function that calls, [2] name of value to test
#'
#' @return TRUE, if no error has occured, stops with error otherwise
test_tweak_vec <- function(tweak.vec,
                           output.info){
  error.message <- paste0("In ", output.info[1], ": ", "'", output.info[2], "'")
  if(!is.numeric(tweak.vec)){
    error.message <- paste0(error.message, " is not numeric")
    stop(error.message, call. = FALSE)
  }
  if(any(is.na(tweak.vec))){
    error.message <- paste0(error.message, " includes NA")
    stop(error.message, call. = FALSE)
  }
  return(TRUE)
}

#' test_numeric
#' (All test functions are used for safety checks in the DTD package)
#' @param test.value value to be tested for integer
#' @param output.info [1] function that calls, [2] name of value to test
#' @param min min is tested '<'
#' @param max max is tested '>'
#'
#' @return TRUE, if no error has occured, stops with error otherwise
test_numeric <- function(test.value,
                         output.info,
                         min,
                         max){

  error.message <- paste0("In ", output.info[1], ": ", "'", output.info[2], "'")

  if(!is.numeric(test.value) || length(test.value) != 1){
    error.message <- paste0(error.message, " is not a single integer")
    stop(error.message, call. = FALSE)
  }
  if(test.value < min){
    error.message <- paste0(error.message, " is below minimal value")
    stop(error.message, call. = FALSE)
  }
  if(test.value > max){
    error.message <- paste0(error.message, " is above maximal value")
    stop(error.message, call. = FALSE)
  }
  return(TRUE)
}

#' test_logical
#' (All test functions are used for safety checks in the DTD package)
#' @param test.value value to be tested for integer
#' @param output.info [1] function that calls, [2] name of value to test
#'
#' @return TRUE, or it throws an error
test_logical <- function(test.value,
                         output.info){
  error.message <- paste0("In ", output.info[1], ": ", "'", output.info[2], "'")

  if(any(!is.logical(test.value)) || length(test.value) != 1){
    error.message <- paste0(error.message, " must be a single value, either 'TRUE' or 'FALSE' (not 1 or 0)")
    stop(error.message, call. = FALSE)
  }
  return(TRUE)
}


#' test_c_type
#' (All test functions are used for safety checks in the DTD package)
#' @param test.value value to be tested for integer
#' @param output.info [1] function that calls, [2] name of value to test
#'
#' @return a estimate_c function, or stops with an error
test_c_type <- function(test.value,
                        output.info){
  error.message <- paste0("In ", output.info[1], ": ", "'", output.info[2], "'")
  if(length(test.value) != 1){
    error.message <- paste0(error.message, " provide a single string to 'estimate.c.type'")
    stop(error.message, call. = FALSE)
  }

  if(!test.value %in% c("non_negative", "direct")){
    error.message <- paste0(error.message, " does not match 'non_negative' or 'direct.")
    stop(error.message, call. = FALSE)
  }
  return(NULL)
}

#' test_string
#' (All test functions are used for safety checks in the DTD package)
#' @param test.value value to be tested for 'as.character' usage
#' @param output.info [1] function that calls, [2] name of value to test
#'
#' @return string
test_string <- function(test.value,
                        output.info){
  useable <- try(as.character(test.value), silent = TRUE)
  if(any(grepl(x = useable, pattern = "Error"))){
    message <- paste0("In ", output.info[1], ": ", "'", output.info[2], "'")
    message <- paste0(message, " can not be used 'as.character'")
    message(message)
    return("")
  }
  else{
    return(test.value)
  }
}

