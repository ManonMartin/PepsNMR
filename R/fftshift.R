fftshift1D2D <- function(x) {
  vec <- F
  if (is.vector(x)) {
    x <- vec2mat(x)
    vec <- T
  }
  m <- dim(x)[2]
  p <- ceiling(m/2)
  new_index <- c((p+1):m, 1:p)
  y <- x[, new_index, drop=vec]
}
