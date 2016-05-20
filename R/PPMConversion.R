#' @export PPMConversion
PPMConversion <- function(RawSpect_data, RawSpect_info,
                          shiftHandling=c("cut", "zerofilling","NAfilling", "circular"), thres = 300) {
  begin_info <- beginTreatment("PPMConversion", RawSpect_data, RawSpect_info)
  RawSpect_data <- begin_info[["Signal_data"]]
  RawSpect_info <- begin_info[["Signal_info"]]
  shiftHandling = match.arg(shiftHandling)

  
  findTMSPpeak <- function(ft, thres. = thres) {
    ft = Re(ft)
    
    res = abs(ft)/cumsum(ft)
    
    N = length(ft)
    seuil = thres.*median(res)
    v = which(res[-c(1:3000)] > seuil)[1]#remove first points
    if (is.na(v)){
      warning("No peak found, need to lower the threshold.")
    }
    v = v + 3000
    
    d = which.max(ft[v:(v+N*0.01)]) # recherche dans les 100 points suivants du max
    new.peak = v+d # pic TMSP 
    
    if (names(which.max(ft[v:(v+N*0.01)])) != names(which.max(ft[v:(v+N*0.03)]))){
      warning("the TMSP peak might be located further away, increase the threshold to check.")
    }
    
    return(new.peak)
  }
  
  n <- nrow(RawSpect_data)
  m <- ncol(RawSpect_data)
  # The Sweep Width has to be the same since the column names are the same
  SW <- RawSpect_info[1,"SW"] # Sweep Width in ppm (semi frequency scale in ppm)
  ppmInterval <- SW / m # FIXME divide by two ??
  TMSPpeaks <- apply(RawSpect_data, 1, findTMSPpeak)
  maxpeak <- max(TMSPpeaks)
  minpeak <- min(TMSPpeaks)
  if (shiftHandling == "zerofilling" || shiftHandling == "NAfilling") {
    fill <- NA
    if (shiftHandling == "zerofilling") {
      fill <- 0
    }
    start <- 1 - maxpeak
    end <- m - minpeak
    ppmScale <- (start:end) * ppmInterval
    Spectrum_data <- matrix(fill, nrow=n, ncol=end-start+1,
                            dimnames=list(rownames(RawSpect_data), ppmScale))
    for (i in 1:n) {
      shift <- (1-TMSPpeaks[i]) - start
      Spectrum_data[i,(1+shift):(m+shift)] <- RawSpect_data[i,]
    }
  } else if (shiftHandling == "cut") {
    start <- 1 - minpeak
    end <- m - maxpeak
    ppmScale <- (start:end) * ppmInterval
    Spectrum_data <- matrix(nrow=n, ncol=end-start+1,
                            dimnames=list(rownames(RawSpect_data), ppmScale))
    for (i in 1:n) {
      Spectrum_data[i,] <- RawSpect_data[i,(1+(TMSPpeaks[i]-minpeak)):(m-(maxpeak-TMSPpeaks[i]))]
    }
  } else { # circular
    start <- 1 - maxpeak
    end <- m - maxpeak
    ppmScale <- (start:end) * ppmInterval
    Spectrum_data <- matrix(nrow=n, ncol=end-start+1,
                            dimnames=list(rownames(RawSpect_data), ppmScale))
    for (i in 1:n) {
      shift <- (maxpeak-TMSPpeaks[i])
      Spectrum_data[i,(1+shift):m] <- RawSpect_data[i,1:(m-shift)]
      if (shift > 0) {
        Spectrum_data[i,1:shift] <- RawSpect_data[i,(m-shift+1):m]
      }
    }
  }
  return(endTreatment("PPMConversion", begin_info, Spectrum_data))
}
