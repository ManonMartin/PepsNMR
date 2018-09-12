#' @export InternalReferencing
InternalReferencing <- function(Spectrum_data, Fid_info, method = c("max", "thres"), 
                          range = c("nearvalue", "all", "window"), ppm.value = 0, 
                          direction = "left", shiftHandling = c("zerofilling", "cut", 
                          "NAfilling", "circular"), c = 2, pc = 0.02, fromto.RC = NULL,
                          ppm.ir = TRUE, rowindex_graph = NULL) {
  
  
  
  # Data initialisation and checks ----------------------------------------------
  
  begin_info <- beginTreatment("InternalReferencing", Spectrum_data, Fid_info)
  Spectrum_data <- begin_info[["Signal_data"]]
  Fid_info <- begin_info[["Signal_info"]]

  
  ######## Check input arguments
  
  range <- match.arg(range)
  shiftHandling <- match.arg(shiftHandling)
  method <- match.arg(method)
  plots <- NULL
  
  checkArg(ppm.ir, c("bool"))
  checkArg(unlist(fromto.RC), c("num"), can.be.null = TRUE)
  checkArg(pc, c("num"))
  checkArg(ppm.value, c("num"))
  checkArg(rowindex_graph, "num", can.be.null = TRUE)
  
  # fromto.RC : if range == "window", 
  # fromto.RC defines the spectral window where to search for the peak
  if (!is.null(fromto.RC)) {
    diff <- diff(unlist(fromto.RC))[1:length(diff(unlist(fromto.RC)))%%2 !=0]
    for (i in 1:length(diff)) {
      if (diff[i] >= 0)  {
        fromto <- c(fromto.RC[[i]][2], fromto.RC[[i]][1])
        fromto.RC[[i]] <- fromto
      }
    }
  }
  
  
  # findTMSPpeak function ----------------------------------------------
  # If method == "tresh", findTMSPpeak will find the position of the first 
  # peak (from left or right) which is higher than a predefined threshold 
  # and is computed as: c*(cumulated_mean/cumulated_sd)
  findTMSPpeak <- function(ft, c = 2, direction = "left") {
    ft <- Re(ft)  # extraction de la partie réelle
    N <- length(ft)
    if (direction == "left") {
      newindex <- rev(1:N)
      ft <- rev(ft)
    }
    thres <- 99999
    i <- 1000  # Start at point 1000 to find the peak
    vect <- ft[1:i]
    
    while (vect[i] <= (c * thres)) {
      cumsd <- stats::sd(vect)
      cummean <- mean(vect)
      thres <- cummean + 3 * cumsd
      i <- i + 1
      vect <- ft[1:i]
    }
    if (direction == "left") {
      v <- newindex[i]
    } else {v <- i}
    
    if (is.na(v))  {
      warning("No peak found, need to lower the threshold.")
      return(NA)
    } else  {
      # recherche dans les 1% de points suivants du max trouve pour etre au sommet du
      # pic
      d <- which.max(ft[v:(v + N * 0.01)])
      new.peak <- v + d - 1  # nouveau pic du TMSP si d > 0
      
      if (names(which.max(ft[v:(v + N * 0.01)])) != names(which.max(ft[v:(v + N * 0.03)])))   {
        # recherche dans les 3% de points suivants du max trouve pour eviter un faux
        # positif
        warning("the TMSP peak might be located further away, increase the threshold to check.")
      }
      return(new.peak)
    }
  }
  
  
  # Define the search zone  ----------------------------------------
  
  n <- nrow(Spectrum_data)
  m <- ncol(Spectrum_data)
  
  # The Sweep Width (SW) has to be the same since the column names are the same
  SW <- Fid_info[1, "SW"]  # Sweep Width in ppm 
  ppmInterval <- SW/(m-1)  # size of a ppm interval
  
  # range: How the search zone is defined ("all", "nearvalue" or "window")
  if (range == "all") {
    
    Data <- Spectrum_data
    
  } else { # range = "nearvalue" or "window"
      # Need to define colindex (column indexes) to apply indexInterval on it
      
    if (range == "nearvalue")  {
        
        fromto.RC <- list(c(-(SW * pc)/2 + ppm.value, (SW * pc)/2 + ppm.value))  # automatic fromto values in ppm
        colindex <- as.numeric(colnames(Spectrum_data))
      
        } else {
          # range == "window"
          # fromto.RC is already user-defined
          if (ppm.ir == TRUE)   {
            colindex <- as.numeric(colnames(Spectrum_data))
          } else   {
            colindex <- 1:m
          }
        }

    # index intervals taking into account the different elements in the list fromto.RC
    Int <- vector("list", length(fromto.RC)) 
    for (i in 1:length(fromto.RC))  {
      Int[[i]] <- indexInterval(colindex, from = fromto.RC[[i]][1], 
                                to = fromto.RC[[i]][2], inclusive = TRUE)
    }
    
    # define Data as the cropped spectrum including the index intervals
    # outside the research zone, the intensities are set to the minimal 
    # intensity of the research zone
    
    if (n > 1){
      Data <- apply(Re(Spectrum_data[,unlist(Int)]),1, function(x) rep(min(x), m))
      Data <- t(Data)
      Data[,unlist(Int)] <- Re(Spectrum_data[,unlist(Int)])
    } else {
      Data <- rep(min(Re(Spectrum_data)) ,m)
      Data[unlist(Int)] <- Re(Spectrum_data[unlist(Int)])
    }
    
  }
  
  
  # Apply the peak location search method ('thres' or 'max') on spectra
  # -----------------------------------------------------------------------
  
  if (method == "thres") {
    TMSPpeaks <- apply(Data, 1, findTMSPpeak, c = c, direction = direction)
  } else { # method == "max
    TMSPpeaks <- apply(Re(Data), 1, which.max)
  }
  
  
  # Shift spectra according to the TMSPpeaks found --------------------------------
  # Depends on the shiftHandling
  
  # TMSPpeaks is a column index
  maxpeak <- max(TMSPpeaks) # max accross spectra
  minpeak <- min(TMSPpeaks) # min accross spectra
  
  
  if (shiftHandling %in% c("zerofilling", "NAfilling",  "cut")) {
    fill <- NA
    if (shiftHandling == "zerofilling")  {
      fill <- 0
    }

    start <-  maxpeak - 1
    end <- minpeak - m
    
    # ppm values of each interval for the whole spectral range of the spectral matrix
    ppmScale <- (start:end) * ppmInterval
    
    # check if ppm.value is in the ppmScale interval
    if(ppm.value < min(ppmScale) | ppm.value > max(ppmScale)) {
    warning("ppm.value = ", ppm.value, " is not in the ppm interval [", 
           round(min(ppmScale),2), ",", round(max(ppmScale),2), "], and is set to its default ppm.value 0")
      ppm.value = 0
      }
    
    # if ppm.value != 0, ppmScale is adapted
    ppmScale <- ppmScale + ppm.value
    
    # create the spectral matrix with realigned spectra
    Spectrum_data_calib <- matrix(fill, nrow = n, ncol =  -(end - start) + 1, 
                            dimnames = list(rownames(Spectrum_data), ppmScale))
    
    # fills in Spectrum_data_calib with shifted spectra
    for (i in 1:n)  {
      shift <- (1 - TMSPpeaks[i]) + start
      Spectrum_data_calib[i, (1 + shift):(m + shift)] <- Spectrum_data[i, ]
    }
    
    if (shiftHandling == "cut")  {
      Spectrum_data_calib = as.matrix(stats::na.omit(t(Spectrum_data_calib)))
      Spectrum_data_calib = t(Spectrum_data_calib)
      base::attr(Spectrum_data_calib, "na.action") <- NULL
    }
  
    
 } else {
    # circular
    start <- 1 - maxpeak
    end <- m - maxpeak
    
    ppmScale <- (start:end) * ppmInterval
    
    # check if ppm.value in is the ppmScale interval
    if(ppm.value < min(ppmScale) | ppm.value > max(ppmScale)) {
      warning("ppm.value = ", ppm.value, " is not in the ppm interval [", 
              round(min(ppmScale),2), ",", round(max(ppmScale),2), "], and is set to its default ppm.value 0")
      ppm.value = 0
    }
    
    # if ppm.value != 0, ppmScale is adapted
    ppmScale <- ppmScale + ppm.value
    
    # create the spectral matrix with realigned spectra
    Spectrum_data_calib <- matrix(nrow=n, ncol=end-start+1,
                            dimnames=list(rownames(Spectrum_data), ppmScale))
    
    # fills in Spectrum_data_calib with shifted spectra
    for (i in 1:n) {
      shift <- (maxpeak-TMSPpeaks[i])
      Spectrum_data_calib[i,(1+shift):m] <- Spectrum_data[i,1:(m-shift)]
      if (shift > 0) {
        Spectrum_data_calib[i,1:shift] <- Spectrum_data[i,(m-shift+1):m]
      }
    }
  }
  
  

  
  # Plot of the spectra (depending on rowindex_graph) ---------------------------------------------------
  
  ppm = xstart = value = xend = Legend = NULL # only for R CMD check
  
  
  # with the search zone for TMSP and the location of the peaks just found
  if (!is.null(rowindex_graph)) {
    
    if (range == "window")  {
      if (ppm.ir == TRUE)   {
        fromto <- fromto.RC
      } else  {
        fromto <- list()
        idcol <- as.numeric(colnames(Spectrum_data))
        for (i in 1:length(fromto.RC)) {
          fromto[[i]] <- as.numeric(colnames(Spectrum_data))[fromto.RC[[i]]]
        }
      }
    } else {
      fromto <- fromto.RC
    }
    
    # TMSPloc in ppm
    TMSPloc <- as.numeric(colnames(Spectrum_data))[TMSPpeaks[rowindex_graph]]
    
    # num plot per window
    num.stacked <- min(6, length(rowindex_graph))
    
    # rectanglar bands of color for the search zone
    rects <- data.frame(xstart = sapply(fromto, function(x) x[[1]]), 
                        xend = sapply(fromto, function(x) x[[2]]), 
                        Legend = "Peak search zone and location")
    
    # vlines for TMSP peak
    addlines <- data.frame(rowname = rownames(Spectrum_data)[rowindex_graph],TMSPloc)
    
    nn <- length(rowindex_graph)
    i <- 1
    j <- 1
    plots <- vector(mode = "list", length = ceiling(nn/num.stacked))
    
    Data <- Spectrum_data[rowindex_graph,]
    while (i <= nn) {
      
      last <- min(i + num.stacked - 1, nn)
      
      melted <- reshape2::melt(Re(Data[i:last, ]), 
                               varnames = c("rowname", "ppm"))
      
      plots[[j]] <- ggplot2::ggplot() + ggplot2::theme_bw() + 
        ggplot2::geom_line(data = melted, 
        ggplot2::aes(x = ppm, y = value)) + 
        ggplot2::geom_rect(data = rects, ggplot2::aes(xmin = xstart, xmax = xend, 
                           ymin = -Inf, ymax = Inf, fill = Legend), alpha = 0.4) + 
        ggplot2::facet_grid(rowname ~ ., scales = "free_y") + 
        ggplot2::theme(legend.position = "none") + 
        ggplot2::geom_vline(data = addlines, ggplot2::aes(xintercept = TMSPloc), 
                            color = "red", show.legend = TRUE) + 
        ggplot2::ggtitle("Peak search zone and location") + 
        ggplot2::theme(legend.position = "top", legend.text = ggplot2::element_text())
      
      
      
      if ((melted[1, "ppm"] - melted[(dim(melted)[1]), "ppm"]) > 0) {
        plots[[j]] <- plots[[j]] + ggplot2::scale_x_reverse()
      }
      
      i <- last + 1
      j <- j + 1
    }
    
    plots
  }
  
  
  # Return the results ----------------------------------------------
  Spectrum_data <- endTreatment("InternalReferencing", begin_info, Spectrum_data_calib)
  
  if (is.null(plots)) {
    return(Spectrum_data)
  } else {
    return(list(Spectrum_data = Spectrum_data, plots = plots))
  }
  
}
