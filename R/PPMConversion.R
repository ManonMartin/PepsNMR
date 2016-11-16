#' @export PPMConversion
PPMConversion <- function(RawSpect_data, RawSpect_info,
                          shiftHandling=c("cut", "zerofilling","NAfilling", "circular"), c = 2, 
                          freq = TRUE, fromto.TMSP = NULL) {
  begin_info <- beginTreatment("PPMConversion", RawSpect_data, RawSpect_info)
  RawSpect_data <- begin_info[["Signal_data"]]
  RawSpect_info <- begin_info[["Signal_info"]]

### CHECK INPUT ARGUMENTS  
  shiftHandling = match.arg(shiftHandling)
  checkArg(freq, c("bool"))
  checkArg(unlist(fromto.TMSP), c("num"), can.be.null = TRUE)
  
  # fromto.TMSP
  diff = diff(unlist(fromto.TMSP))[1:length(diff(unlist(fromto.TMSP))) %% 2 != 0]
  for (i in 1:length(diff)) {
    if (diff[i]<=0) {
      stop(paste("Invalid region removal because from > to"))
    } 
  }
  
  
################################## 
### findTMSPpeak FUNCTION
##################################  
  findTMSPpeak <- function(ft, c=2) {
    ft = Re(ft) # extraction de la partie réelle
    N = length(ft)
    thres = 99999
    i= 1000 # Start at point 1000 to find the peak
    vect = ft[1:i]
    while (vect[i] <= (c*thres)) {
      cumsd = stats::sd(vect)
      cummean = mean(vect)
      thres = cummean + 3*cumsd
      i=i+1
      vect = ft[1:i]
    }
    v = i
    if (is.na(v)){
      warning("No peak found, need to lower the threshold.")
      return(NA)
    } else{
      d = which.max(ft[v:(v+N*0.01)]) # recherche dans les 1% de points suivants du max trouve pour etre au sommet du pic
      new.peak = v+d-1 # nouveau pic du TMSP si d > 0
      
      if (names(which.max(ft[v:(v+N*0.01)])) != names(which.max(ft[v:(v+N*0.03)]))){ # recherche dans les 3% de points suivants du max trouve pour eviter un faux positif
        warning("the TMSP peak might be located further away, increase the threshold to check.")
      }
      return(new.peak)
    }
  }
  
  
################################## 
### APPLY findTMSPpeak ON SPECTRA
################################## 
  
  n <- nrow(RawSpect_data)
  m <- ncol(RawSpect_data)
  # The Sweep Width has to be the same since the column names are the same
  SW <- RawSpect_info[1,"SW"] # Sweep Width in ppm (semi frequency scale in ppm)
  ppmInterval <- SW / m # FIXME divide by two ??
  
  if (is.null(fromto.TMSP)) {
    Data = RawSpect_data
  }else{
    
  # if freq == TRUE, then fromto is in the colnames values, else, in the column index
    if (freq == TRUE) {
      colindex = as.numeric(colnames(RawSpect_data))
    }else{ colindex = 1:m}
    
    Int = vector("list", length(fromto.TMSP))
    for (i in 1:length(fromto.TMSP)) {
      Int[[i]] <- indexInterval(colindex, 
                                from = fromto.TMSP[[i]][1], 
                                to = fromto.TMSP[[i]][2], inclusive=TRUE)
    }
    
    vector = rep(0, m)
    vector[unlist(Int)]=1
    if (n>1) {
      Data = sweep(RawSpect_data, MARGIN=2,  FUN = "*", vector) # Cropped_Spectrum
    } else {Data = RawSpect_data*vector} # Cropped_Spectrum
  }
  
  
  TMSPpeaks <- apply(Data, 1, findTMSPpeak)
  maxpeak <- max(TMSPpeaks)
  minpeak <- min(TMSPpeaks)
  
  
################################## 
### SHIFT SPECTRA
##################################   
  
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
  
  
################################## 
### RETURN RESULTS
################################## 
  
  return(endTreatment("PPMConversion", begin_info, Spectrum_data))
}
