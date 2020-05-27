#' @title BotrytisConverter
#'
#' @description Converts the percentage scores (0-100%) of botrytis into character (0-9) scores using the Seges Septoria conversion scale
#'
#' @param perc
#'
#' return value1
#'
#' @examples
#' newvalue = conversion_func_threshold(0.5)
#'
#' @export

conversion_func_threshold <- function(perc){
  # A function that takes a percentage botrytis score
  # and translate to the character in the dictionary that
  # is closest to. If the percentage is equally close
  # to two values, chose the highest number (more susceptibility)
  h <- hash()
  h[["0"]] <- 0
  h[["0.01"]] <- 1
  h[["0.1"]] <- 2
  h[["0.5"]] <- 3
  h[["1"]] <- 4
  h[["5"]] <- 5
  h[["10"]] <- 6
  h[["25"]] <- 7
  h[["75"]] <- 8
  h[["100"]] <- 9

  if (is.na(perc)==T){
    return(perc)
  }

  for (i in seq(1:length(keys(h)))){
    candidates=keys(h)[which(as.numeric(keys(h))<=perc)]
    best_translation=max(as.numeric(candidates))
  }

  value1=h[[as.character(best_translation)]]

  return(value1)
}

