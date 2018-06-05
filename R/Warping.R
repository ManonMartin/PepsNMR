#' @export Warping
Warping <- function(Spectrum_data, normalization.type = c("median", "mean", 
                    "firstquartile", "peak", "none"), fromto.normW = c(3.05, 4.05), 
                    reference.choice = c("fixed", "before", "after", "manual"), reference = 1, 
                    optim.crit = c("RMS", "WCC"), ptw.wp = F, K = 3, L = 40, lambda.smooth = 0, 
                    deg = 3, lambda.bspline = 0.01, kappa = 1e-04, max_it_Bspline = 10, 
                    returnReference = FALSE, returnWarpFunc = FALSE) {
  
  # Mean square difference function definition --------------------------------------
  
  meanSqrDiff <- function(m, row){
    # for the row ref, x - m[ref,] is 0 to get, the mean, we divide by
    # nrow(data)-1 because the row ref is ignored we could actually just do
    # the sum and not the mean, both choices are good
    return(sum(apply(m, 1, function(x) sum((x - m[row, ])^2, na.rm = T)))/(nrow(m) - 
      1))
  }
  
  # Data initialisation and checks ----------------------------------------------

  begin_info <- beginTreatment("Warping", Spectrum_data, force.real = T)
  Spectrum_data <- begin_info[["Signal_data"]]
  normalization.type <- match.arg(normalization.type)
  reference.choice <- match.arg(reference.choice)
  optim.crit <- match.arg(optim.crit)
  
  checkArg(K, c("int", "pos"))
  checkArg(L, c("int", "pos0"))
  checkArg(lambda.smooth, c("num", "pos0"))
  checkArg(deg, c("int", "pos"))
  checkArg(lambda.bspline, c("num", "pos0"))
  checkArg(kappa, c("num", "pos0"))
  checkArg(max_it_Bspline, c("int", "pos"))
  
  if (reference.choice == "fixed" & !(reference[1] %in% row.names(Spectrum_data))) {
    checkArg(reference, c("int", "pos"))
    reference <- row.names(Spectrum_data)[reference]
  }
  
  if (reference.choice == "manual" & (length(reference) != dim(Spectrum_data)[2] | !is.numeric(reference))) {
  stop("reference is misspecified with reference.choice == manual")
  }
  
  if (L > 0 && L <= deg) {
    stop("L should be greater than deg because with 1 interval, there is already deg+1 Bsplines.")
  }
 
  
  # adding the manual reference to the spectral matrix 
  if (reference.choice == "manual") {
    colnam <- colnames(Spectrum_data)
    rownam <- rownames(Spectrum_data)
    Spectrum_data <- rbind(reference, Spectrum_data)
    Spectrum_data <- as.matrix(Spectrum_data)
    colnames(Spectrum_data) <- colnam
    rownames(Spectrum_data) <- c("manual_ref", rownam)
  }
  
  n <- nrow(Spectrum_data)
  m <- ncol(Spectrum_data)
  
  # Data pre-normalization ----------------------------------------------
  
  if (normalization.type != "none") {
    norm.res <- Normalization(Spectrum_data, type.norm = normalization.type, 
                                   fromto.norm = fromto.normW, returnFactor = TRUE)
    Spectrum_data <- norm.res[["Spectrum_data"]]
    normFactor <- norm.res[["factor"]]
  }
  
  # Warping -----------------------------------------------------
  
  rnames <- rownames(Spectrum_data)
  if (n > 1) {
    # (pool of potential) candidate(s) as reference spectrum
    if (reference.choice == "fixed") {
      pool <- reference
    } else if (reference.choice == "before") {
      argmin <- which.min(sapply(rnames, function(i) meanSqrDiff(Spectrum_data, 
        i)))
      pool <- names(argmin)
    } else if (reference.choice == "after") {
      pool <- rnames
    } else {pool = "manual_ref"} # dummy value that works
    
    best.meanSqrDiff <- NULL
    best.Warped_data <- NULL
    decreasing <- FALSE
    

    
    # Warping for each potential reference specturm
      for (reference in pool) {

        ref <- Spectrum_data[reference, ]  # reference spectrum 
        
        samp_rnames <- rnames[rnames != reference]
        sample <- Spectrum_data[samp_rnames, , drop = F]  # spectra to be warped
        beta <- rep(0, K + 1)  # starting coefficients
        cur.Warped_data <- Spectrum_data  # initialize the matrix of warped spectra
        warp.func <- Spectrum_data  # initialize the matrix of estimated warping functions
        
        # if (reference.choice %in% c("fixed", "before", "after")) { 
        warp.func[reference, ] <- 0
        # }
        
        if (K >= 1) {
          # This is beta_1 which should be approximately 1 because if they are
          # perfectly aligned w(v) = v = 0*1 + 1*v + 0*v^2 + 0*v^3
          beta[2] <- 1
        }
        if (ptw.wp)  {
          
          # A. Parametric time warping with polynomials
          ptw.output <- ptw::ptw(ref, sample, optim.crit = optim.crit, 
            init.coef = beta, smooth.param = lambda.smooth)
          cur.Warped_data[samp_rnames, ] <- ptw.output$warped.sample
          # We transpose because 'diff' is along the columns
          w <- t(ptw.output$warp.fun)
          warp.func[samp_rnames, ] <- w
        } else  {
          
          # or B. Semi-parametric time warping with polynomials AND B splines
          if (optim.crit == "WCC")  {
            stop("WCC is only implemented in ptw, set ptw.wp=T to use WCC.")
          }
          for (samp_rname in samp_rnames)  {
            sw.output <- SingleWarp(ref = ref, sample = sample[samp_rname, 
            ], beta = beta, L = L, lambda.smooth = lambda.smooth, 
            deg = deg, lambda.bspline = lambda.bspline, kappa = kappa, 
            max_it_Bspline = max_it_Bspline)
            warped <- sw.output$warped
            w <- sw.output$w
            cur.Warped_data[samp_rname, ] <- warped
            warp.func[samp_rname, ] <- w
          }
        }
        cur.meanSqrDiff <- meanSqrDiff(cur.Warped_data, reference)
        
        # Update the results if current is best
        if (is.null(best.meanSqrDiff) || cur.meanSqrDiff < best.meanSqrDiff) {
          best.meanSqrDiff <- cur.meanSqrDiff
          best.Warped_data <- cur.Warped_data
          decreasing <- (FALSE %in% (diff(w) > 0))
        }
      }

    
    # check if warping function is decreasing
    if (decreasing) {
      warning("The warping function is not increasing for the sample", 
        samp_rname, ".")
    }
  } else {
    best.Warped_data <- Spectrum_data
  }
  
  
  
  # Data finalisation ----------------------------------------------
  best.Warped_data <- best.Warped_data*normFactor # denormalize the spectra
  
  if (reference.choice=="manual") {
    best.Warped_data <- best.Warped_data[-1,] # remove the manual reference spectrum
  }
  
  Spectrum_data <- endTreatment("Warping", begin_info, best.Warped_data)
  if (returnReference) {
    if (returnWarpFunc)  {
      return(list(Spectrum_data = Spectrum_data, Reference = reference, 
        Warp.func = warp.func))
    } else {
      return(list(Spectrum_data = Spectrum_data, Reference = reference))
    }
  } else {
    if (returnWarpFunc) {
      return(list(Spectrum_data = Spectrum_data, Warp.func = warp.func))
    } else  {
      return(Spectrum_data)
    }
  }
  
}
