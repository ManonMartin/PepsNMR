# Shift zero-frequency component to center of spectrum
# When applying FT for discrete data, the multiplication 
# shifts the spectrum of the signal by half the sampling frequency
# ==> the zero frequency that was at 0 is now at half width

fftshift1D2D <- function(x) {
  vec <- FALSE
  if (is.vector(x)) {
    x <- vec2mat(x)
    vec <- TRUE
  }
  m <- dim(x)[2]
  p <- ceiling(m/2)
  new_index <- c((p + 1):m, 1:p)
  y <- x[, new_index, drop = vec]
}
