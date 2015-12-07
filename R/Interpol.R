Interpol <- function(t, y) {
  m <- length(y)
  # t <= m-1
  # because if t > m-1, y[ti+1] will be NA when we compute g
  valid <- 1 <= t & t <= m-1 # FIXME it was '<' in Bubble v2
  s <- (1:m)[valid]
  ti <- floor(t[s])
  tr <- t[s] - ti
  g <- y[ti + 1] - y[ti]
  f <- y[ti] + tr * g
  list(f=f, s=s, g=g)
}