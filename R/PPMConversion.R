#' @export PPMConversion
PPMConversion <- function(RawSpect_data, RawSpect_info,
                          shiftHandling=c("cut", "zerofilling","NAfilling", "circular"), 
                          from.ppmc = 7400, to.ppmc = 9400) {
  begin_info <- beginTreatment("PPMConversion", RawSpect_data, RawSpect_info)
  RawSpect_data <- begin_info[["Signal_data"]]
  RawSpect_info <- begin_info[["Signal_info"]]
  shiftHandling = match.arg(shiftHandling)
  checkArg(from.ppmc, "int", "pos")
  checkArg(to.ppmc, "int", "pos")
  
  findTMSPpeak <- function(ft) {
    # The range is from 7400 to 9400 when we have 2^15 points
    N <- length(ft)
    M=2^15
    FROM <- floor((from * N) / M)
    TO <- ceiling((to * N) / M)
    searchZone <- Re(ft[FROM:TO])
    peakInZone <- which.max(searchZone)
    # peakInZone is a numeric with a name
    # The name is the time value of the max
    # the numerical value is the index of the max
    return(FROM-1 + peakInZone)
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
