PolyWarp <- function(ref, sample, beta) {
  
  ncoef <- length(beta)
  
  m <- max(length(ref), length(sample))
  
  # Divide by m to avoid too large numbers (when we take the cube for example)
  time <- (1:m)/m
  
  # We need to adapt beta with our change
  a <- beta * m^(0:(ncoef - 1))  # polynomial coefficients
  
  B <- matrix(time, nrow = m, ncol = ncoef)  # B = t^k
  B <- t(apply(B, 1, cumprod))/B
  
  rms_old <- 0
  da <- 0  # update the coefficients
  for (it in 1:40) {
    # I do it here to avoir incrementing a at the last iteration
    a <- a + da
    # Compute warping function, w
    w <- B %*% a
    # Warp template x by linear interpolation, to give z = x(w(x)) and its derivative
    # g
    interp.out <- Interpol(w, sample)
    z <- interp.out$f
    sel <- interp.out$s
    g <- interp.out$g
    
    # Compute residuals and check convergence
    r <- ref[sel] - z
    rms <- sqrt(mean(r * r, na.rm = T))
    drms <- abs((rms - rms_old)/(rms + 1e-10))
    if (drms < 1e-06) {
      break
    }
    rms_old <- rms
    
    # Improve coeffcients with linear regression
    G <- kronecker(matrix(1, nrow = 1, ncol = ncoef), g)
    Q <- G * B[sel, ]
    da <- qr.solve(Q, r)
  }
  
  ## back-transform into beta coefficients
  a <- a/m^(0:(ncoef - 1))
  
  # ptw::interpol(w, sample)
  list(w = w, sel = sel, beta = a, warped = z)
}
