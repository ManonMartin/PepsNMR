#' @export WindowSelection
WindowSelection <- function(Spectrum_data, from.ws = 10, to.ws = 0.2) {
  
  
  # Data initialisation and checks ----------------------------------------------
  
  begin_info <- beginTreatment("WindowSelection", Spectrum_data)
  Spectrum_data <- begin_info[["Signal_data"]]
  checkArg(from.ws, "num", can.be.null = TRUE)
  checkArg(to.ws, "num", can.be.null = TRUE)
  
  m <- ncol(Spectrum_data)
  
  # largestWindowWithoutNA function  ----------------------------------------------
  # Used to find the global largest window without NA
  largestWindowWithoutNA <- function(data, from_index = NULL, to_index = NULL) {
    largestFromOrToWithoutNA <- function(fromorto_index, fromorto, delta)  {
      if (is.na(data[fromorto_index]))  {
        stop("I should return an empty matrix because there is NA at" + fromorto + "ppm")
      }
      toorfrom_index <- fromorto_index
      while (toorfrom_index > 0 && toorfrom_index <= m && !is.na(data[toorfrom_index])) {
        toorfrom_index <- toorfrom_index + delta
      }
      toorfrom_index <- toorfrom_index - delta
      return(toorfrom_index)
    }
    if (!is.null(from_index) && !is.null(to_index))  {
      if (is.na(sum(data[from_index:to_index]))) {
        warning("There is NA in the selected window")
      }
    } else if (!is.null(from_index)) {
      to_index <- largestFromOrToWithoutNA(from_index, from.ws, 1)
    } else if (!is.null(to_index)) {
      from_index <- largestFromOrToWithoutNA(to_index, to.ws, -1)
    } else  {
      # Largest interval without NA
      maxLen <- 0
      curLen <- 0
      for (i in 1:m)   {
        if (is.na(data[i]))    {
          curLen <- 0
        } else  {
          curLen <- curLen + 1
          if (curLen > maxLen)
          {
          maxLen <- curLen
          from_index <- i - curLen + 1
          to_index <- i
          }
        }
      }
    }
    return(from_index:to_index)
  }
  
  
  # Window Selection  ----------------------------------------------
  ppm <- as.numeric(colnames(Spectrum_data))

  from_index <- NULL
  to_index <- NULL
  if (!is.null(from.ws)) {
    from_index <- binarySearch(ppm, from.ws, TRUE)
  }
  if (!is.null(to.ws)) {
    to_index <- binarySearch(ppm, to.ws, FALSE)
  }
  interval <- largestWindowWithoutNA(colSums(Spectrum_data), from_index, to_index)
  
  # Data finalisation ----------------------------------------------
  return(endTreatment("WindowSelection", begin_info, Spectrum_data[, interval, 
    drop = FALSE]))
}
