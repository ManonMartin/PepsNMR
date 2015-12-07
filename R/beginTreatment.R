beginTreatment <- function (name, Signal_data=NULL, Signal_info=NULL, force.real=FALSE) {
  cat("Begin", name, "\n")
  vec <- is.vector(Signal_data)
  if (vec) {
    Signal_data <- vec2mat(Signal_data)
  }
  if (is.vector(Signal_info)) {
    Signal_info <- vec2mat(Signal_info)
  }
  if (!is.null(Signal_data)) {
    if (!is.matrix(Signal_data)) {
      stop("Signal_data is not a matrix.")
    }
    if (!is.complex(Signal_data) && !is.numeric(Signal_data)) {
      stop("Signal_data contains non-numerical values.")
    }
  }
  if (!is.null(Signal_info) && !is.matrix(Signal_info)) {
    stop("Signal_info is not a matrix.")
  }
  Original_data <- Signal_data
  if (force.real) {
    if (is.complex(Signal_data)) {
      Signal_data <- Re(Signal_data)
    } else {
      # The signal is numeric
      # Im(Signal_data) is zero anyway so let's
      # avoid using complex(real=...,imaginary=0)
      # which would give a complex signal in endTreatment()
      force.real <- FALSE
    }
  }
  return(list(start=proc.time(),
              vec=vec,
              force.real=force.real,
              Original_data=Original_data,
              Signal_data=Signal_data,
              Signal_info=Signal_info))
}
