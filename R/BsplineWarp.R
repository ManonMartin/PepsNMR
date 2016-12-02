BsplineWarp <- function(ref, sample, m, t, w0, L, deg, lambda.bspline, 
              kappa, max_it) {
  
  # w0: # warping function from polynomial warping
  
  
  # B-spline basis for warping function
  B <- Bspline(t, x0 = 0, x1 = m, ndx = L - deg, deg = deg)
  nb <- ncol(B)
  
  # Penalty matrices
  P <- 0
  if (lambda.bspline > 0 || kappa > 0) {
    I <- diag(1, nb)
    P <- kappa * I
    if (lambda.bspline > 0) {
      D <- diff(I, 3)
      P <- P + lambda.bspline * t(D) %*% D
    }
  }
  
  # Initialize;
  a <- matrix(0, nrow = nb, ncol = 1)  # B-Spline coefficients
  for (it in 1:max_it) {
    # Interpolate
    w <- w0 + B %*% a
    interp.out <- Interpol(w, sample)
    z <- interp.out$f
    sel <- interp.out$s
    g <- interp.out$g
    
    # Improve coefficients
    r <- ref[sel] - z
    # make matrix with nb identical columns 'g'
    G <- kronecker(matrix(1, nrow = 1, ncol = nb), g)
    Q <- G * B[sel, ]
    a <- solve(t(Q) %*% Q + P, t(Q) %*% (r + Q %*% a))
  }
  return(list(w = w, sel = sel, alpha = a, warped = z))
}
