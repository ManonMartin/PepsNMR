Bspline <- function (x, x0, x1, ndx, deg=3) {
  if (min(x) < x0 | max(x) > x1) {
    stop('Some elements of x out of bounds !!')
  }
  dx = (x1 - x0) / ndx
  t = x0 + dx * ((-deg):(ndx-1))
  Tt = matrix(1, nrow=length(x), ncol=1) %*% t
  X = x %*% matrix(1, nrow=1, ncol=length(t))
  D = (X - Tt) / dx
  B = t(diff(rbind(rep.int(0, length(x)), t(D <= 1))))
  r = c(2:length(t), 1)
  for (k in 1:deg) {
    B = (D * B + (k + 1 - D) * B[, r]) / k
  } #FIXME sparse ?
  return(B)
}