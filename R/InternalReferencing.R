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

  
  # Check input arguments
  range <- match.arg(range)
  shiftHandling <- match.arg(shiftHandling)
  method <- match.arg(method)
  plots <- NULL
  
  
  checkArg(ppm.ir, c("bool"))
  checkArg(unlist(fromto.RC), c("num"), can.be.null = TRUE)
  checkArg(pc, c("num"))
  checkArg(ppm.value, c("num"))
  checkArg(rowindex_graph, "num", can.be.null = TRUE)
  
  # fromto.RC
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
  
  
  # Apply the method ('thres' or 'max') on spectra
  # ----------------------------------------------
  
  n <- nrow(Spectrum_data)
  m <- ncol(Spectrum_data)
  
  # The Sweep Width has to be the same since the column names are the same
  SW <- Fid_info[1, "SW"]  # Sweep Width in ppm (semi frequency scale in ppm)
  ppmInterval <- SW/(m-1)  
  
  if (range == "all") {
    Data <- Spectrum_data
  } else { # range = "nearvalue" or "window"
      
      if (range == "nearvalue")  {
        fromto.RC <- list(c(-(SW * pc)/2 + ppm.value, (SW * pc)/2 + ppm.value))  # automatic fromto values in ppm
        colindex <- as.numeric(colnames(Spectrum_data))
      } else {
          if (ppm.ir == TRUE)   {
            colindex <- as.numeric(colnames(Spectrum_data))
          } else   {
            colindex <- 1:m
          }
        }

    
    Int <- vector("list", length(fromto.RC))
    for (i in 1:length(fromto.RC))  {
      Int[[i]] <- indexInterval(colindex, from = fromto.RC[[i]][1], 
                                to = fromto.RC[[i]][2], inclusive = TRUE)
    }
    
    vector <- rep(0, m)
    vector[unlist(Int)] <- 1
    if (n > 1)  {
      Data <- sweep(Spectrum_data, MARGIN = 2, FUN = "*", vector)  # Cropped_Spectrum
    } else  {
      Data <- Spectrum_data * vector
    }  # Cropped_Spectrum
  }
  
  
  if (method == "thres") {
    TMSPpeaks <- apply(Data, 1, findTMSPpeak, c = c, direction = direction)
  } else {
    TMSPpeaks <- apply(abs(Re(Data)), 1, which.max)
  }
  
  # TMSPpeaks is an column index
  maxpeak <- max(TMSPpeaks)
  minpeak <- min(TMSPpeaks)
  
  
  
  # Shift spectra according to the TMSPpeaks found --------------------------------
  # Depends on the shiftHandling
  
  if (shiftHandling %in% c("zerofilling", "NAfilling",  "cut")) {
    fill <- NA
    if (shiftHandling == "zerofilling")  {
      fill <- 0
    }

    start <-  maxpeak - 1
    end <- minpeak - m
    
    ppmScale <- (start:end) * ppmInterval
    
    # check if ppm.value is in the ppmScale interval
    if(ppm.value < min(ppmScale) | ppm.value > max(ppmScale)) {
    warning("ppm.value = ", ppm.value, " is not in the ppm interval [", 
           round(min(ppmScale),2), ",", round(max(ppmScale),2), "], and is set to its default ppm.value 0")
      ppm.value = 0
      }
    
    ppmScale <- ppmScale + ppm.value
    
    Spectrum_data_calib <- matrix(fill, nrow = n, ncol =  -(end - start) + 1, 
                            dimnames = list(rownames(Spectrum_data), ppmScale))
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
    ppmScale <- ppmScale + ppm.value
    
    Spectrum_data_calib <- matrix(nrow=n, ncol=end-start+1,
                            dimnames=list(rownames(Spectrum_data), ppmScale))
    for (i in 1:n) {
      shift <- (maxpeak-TMSPpeaks[i])
      Spectrum_data_calib[i,(1+shift):m] <- Spectrum_data[i,1:(m-shift)]
      if (shift > 0) {
        Spectrum_data_calib[i,1:shift] <- Spectrum_data[i,(m-shift+1):m]
      }
    }
  }
  
  

  
  # Plot of the spectra ---------------------------------------------------
  
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
    num.stacked <- 6
    
    # rectanglar bands of color for the search zone
    rects <- data.frame(xstart = sapply(fromto, function(x) x[[1]]), 
                        xend = sapply(fromto, function(x) x[[2]]), 
                        Legend = "TMSP search zone and location")
    
    # vlines for TMSP peak
    addlines <- data.frame(rowname = rownames(Spectrum_data)[rowindex_graph],TMSPloc)
    
    nn <- length(rowindex_graph)
    i <- 1
    j <- 1
    plots <- vector(mode = "list", length = ceiling(nn/num.stacked))
    
    while (i <= nn) {
      
      last <- min(i + num.stacked - 1, nn)
      
      melted <- reshape2::melt(Re(Spectrum_data[i:last, ]), 
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
        ggplot2::ggtitle("TMSP peak search zone and location") + 
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
