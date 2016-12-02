SingleWarp <- function(ref, sample, beta, L = 40, lambda.smooth = 0, deg = 3, 
                        lambda.bspline, kappa, max_it_Bspline) {
  # specific parameters: for polynomial warping: beta for B spline
  # warping: m, t, w0, L, deg, lambda.bspline, kappa, max_it
  
  m <- length(ref)
  t <- 1:m
  sample0 <- sample
  sel0 <- 1:m
  # e.g. if you don't warp the extremities, you can set sel0 =
  # 7700:25000;
  t <- t[sel0]
  t <- t - min(t) + 0.5
  ref <- ref[sel0]
  sample <- sample[sel0]
  m <- length(ref)
  # From here, we stop using sel0.  We just ignore t, ref and sample out
  # of sel0.
  
  # POLYNOMIAL WARPING -------------------------------------------------
  
  # smoothing
  if (lambda.smooth > 0) {
    # Start with parametric warp, after heavy smoothing
    smooth.sample <- ptw::difsm(sample, lambda.smooth)
    smooth.ref <- ptw::difsm(ref, lambda.smooth)
  } else {
    smooth.sample <- sample
    smooth.ref <- ref
  }
  
  # Warp with polynomials, 1, v, ..., v^K
  pw.out <- PolyWarp(ref = smooth.ref, sample = smooth.sample, beta = beta)
  w <- pw.out$w  # warping function
  sel <- pw.out$sel  # selected indices in spectrum
  beta <- pw.out$beta  # polynomial coefficients
  
  if (lambda.smooth > 0) {
    # We warp the non-smooth sample
    interp.out <- Interpol(w, sample)
    sel <- interp.out$s
    sample[sel] <- interp.out$f
  } else {
    # It has not been smoothed
    sample[sel] <- pw.out$warped
  }
  
  if (L >= 1) {
    
    # BSPLINE WARPING -------------------------------------------------
    bw.out <- BsplineWarp(ref = ref, sample = sample, m = m, t = t, 
      w0 = w, L = L, deg = deg, lambda.bspline = lambda.bspline, 
      kappa = kappa, max_it = max_it_Bspline)
    w <- bw.out$w
    sel <- bw.out$sel
    alpha <- bw.out$alpha
    sample[sel] <- bw.out$warped
  }
  
  # Change the part of sample we have not ignored
  sample0[sel0] <- sample
  return(list(warped = sample0, w = w))
}
