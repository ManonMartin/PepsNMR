#' @export Apodization
Apodization <- function(Fid_data, Fid_info = NULL, DT = NULL, 
                        type.apod = c("exp","cos2", "blockexp", "blockcos2", 
                        "gauss", "hanning", "hamming"), phase = 0, rectRatio = 1/2, 
                        gaussLB = 1, expLB = 1, plotWindow = F, returnFactor = F) {
  
    # Data initialisation and checks ----------------------------------------------
    begin_info <- beginTreatment("Apodization", Fid_data, Fid_info)
    Fid_data <- begin_info[["Signal_data"]]
    Fid_info <- begin_info[["Signal_info"]]
    # Data check
    type.apod <- match.arg(type.apod)
    checkArg(DT, c("num", "pos"), can.be.null = TRUE)
    checkArg(phase, c("num"))
    
    # Apodization ----------------------------------------------
    DT <- getArg(DT, Fid_info, "DT")  # Dwell Time
    m <- ncol(Fid_data)
    t <- (1:m) * DT  # Time
    rectSize <- ceiling(rectRatio * m)
    gaussLB <- (gaussLB/(sqrt(8 * log(2))))
    # Define the types of apodization:
    switch(type.apod, exp = {
        # exponential
        Factor <- exp(-expLB * t)
    }, cos2 = {
        # cos^2
        c <- cos((1:m) * pi/(2 * m) - phase * pi/2)
        Factor <- c * c
    }, blockexp = {
        # block and exponential
        Factor <- c(rep.int(1, rectSize), rep.int(0, m - rectSize))
        # | rectSize | 1 ___________ | \ 0 \____
        Factor[(rectSize + 1):m] <- exp(-expLB * t[1:(m - rectSize)])
    }, blockcos2 = {
        # block and cos^2
        Factor <- c(rep.int(1, rectSize), rep.int(0, m - rectSize))
        c <- cos((1:(m - rectSize)) * pi/(2 * (m - rectSize)))
        Factor[(rectSize + 1):m] <- c * c
    }, gauss = {
        # gaussian
        Factor <- exp(-(gaussLB * t)^2/2)
        Factor <- Factor/max(Factor)
    }, hanning = {
        # Hanning
        Factor <- 0.5 + 0.5 * cos((1:m) * pi/m - phase * pi)
    }, hamming = {
        # Hamming
        Factor <- 0.54 + 0.46 * cos((1:m) * pi/m - phase * pi)
    })
    if (plotWindow) {
        graphics::plot(1:m, Factor, "l")
        # dev.off() # device independent, it is the responsability of the
        # caller to do it
    }
    # Apply the apodization factor on the spectra
    Fid_data <- sweep(Fid_data, MARGIN = 2, Factor, `*`)
    
    # Data finalisation ----------------------------------------------
    Fid_data <- endTreatment("Apodization", begin_info, Fid_data)
    if (returnFactor) {
        return(list(Fid_data = Fid_data, Factor = Factor))
    } else {
        return(Fid_data)
    }
}
