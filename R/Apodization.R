Apodization <- function(Fid_data, Fid_info=NULL, DT=NULL,
                        type = c("exp", "cos2", "blockexp", "blockcos2", "gauss", "hanning", "hamming"),
                        phase=0, rectRatio=1/2, gaussLB=1, expLB=1, plotWindow=F, returnFactor=F) {
  begin_info <- beginTreatment("Apodization", Fid_data, Fid_info)
  Fid_data <- begin_info[["Signal_data"]]
  Fid_info <- begin_info[["Signal_info"]]
  type <- match.arg(type)
  checkArg(DT, c("num", "pos"), can.be.null=TRUE)
  checkArg(phase, c("num"))
  DT <- getArg(DT, Fid_info, "DT") # Dwell Time
  m <- ncol(Fid_data)
  t <- (1:m) * DT # Time
  rectSize <- ceiling(rectRatio*m)
  gaussLB <- (gaussLB/(sqrt(8*log(2))))
  switch(type,
         "exp" = { # exponential
           factor <- exp(-expLB*t)
         },
         "cos2" = { # cos^2
           c <- cos((1:m)*pi/(2*m) - phase*pi/2)
           factor <- c * c
         },
         "blockexp" = { # block and exponential
           factor <- c(rep.int(1, rectSize), rep.int(0, m-rectSize))
           #   | rectSize |
           # 1 ___________
           #              |
           #               \
           # 0              \____
           factor[(rectSize+1):m] <- exp(-expLB*t[1:(m-rectSize)])
         },
         "blockcos2" = { # block and cos^2
           factor <- c(rep.int(1, rectSize), rep.int(0, m-rectSize))
           c <- cos((1:(m-rectSize))*pi/(2*(m-rectSize)))
           factor[(rectSize+1):m] <- c * c
         },
         "gauss" = { # gaussian
           factor <- exp(-(gaussLB*t)^2/2)
           factor <- factor / max(factor)
         },
         "hanning" = { # Hanning
           factor <- 0.5 + 0.5*cos((1:m)*pi/m - phase*pi)
         },
         "hamming" = { # Hamming
           factor <- 0.54 + 0.46*cos((1:m)*pi/m - phase*pi)
         }
  )
  if (plotWindow) {
    plot(1:m, factor, "l")
    # dev.off() # device independent, it is the responsability of the caller to do it
  }
  Fid_data <- sweep(Fid_data, MARGIN=2, factor, `*`)
  Fid_data <- endTreatment("Apodization", begin_info, Fid_data)
  if (returnFactor) {
    return(list(Fid_data=Fid_data, factor=factor))
  } else {
    return(Fid_data)
  }
}
