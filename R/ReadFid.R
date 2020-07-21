#' @export ReadFid
ReadFid <- function(path) {
  
  # Read 1D FID using Bruker XWinNMR and TopSpin format.  It is inspired of the
  # matNMR matlab library which deals with 2D FID and also other formats
  # Read also the parameters in the acqus file
  
  paramFile <- file.path(path, "acqus")
  # BYTEORDA: 0 -> Little Endian 1 -> Big Endian
  params <- readParams(paramFile, c("TD", "BYTORDA", "DIGMOD", "DECIM", "DSPFVS", 
                                    "SW_h", "SW", "O1"))
  
  if (params[["DSPFVS"]] >= 20) {
    # The group delay first order phase correction is given directly from version 20
    grpdly <- readParams(paramFile, c("GRPDLY"))
    params[["GRPDLY"]] <- grpdly[["GRPDLY"]]
  }
  TD <- params[["TD"]]
  endianness <- if (params$BYTORDA) 
    "big" else "little"
  if (TD%%2 != 0) {
    stop(paste("Only even numbers are allowed for size in TD because it is complex 
               data with the real and imaginary part for each element.", 
               "The TD value is in the", paramFile, "file"))
  }
  
  # Interpret params Dwell Time, time between 2 data points in the FID
  params[["DT"]] <- 1/(2 * params[["SW_h"]])
  
  # determine data type (assume int32 unless DTYPA is 2)
  isfloat = FALSE     # default value
  if (params[["DTYPA"]] == 2){
    isfloat = TRUE
  } 
  
  # Read fid
  fidFile <- file.path(path, "fid")
  
  if (isfloat==FALSE){
    fidOnDisk <- readBin(fidFile, what = "int", n = TD, size = 4L, endian = endianness)
  }else {
    fidOnDisk <- readBin(fidFile, what = "double", n = TD, endian = endianness)
  }
  
  # Real size that is on disk (it should be equal to TD2, except for TopSpin/Bruker
  # (which is our case) according to matNMR as just discussed
  TDOnDisk <- length(fidOnDisk)
  if (TDOnDisk < TD) {
    warning("Size is smaller than expected, the rest is filled with zero so the size is the same for every fid")
    fidGoodSize <- sapply(vector("list", length = TD), function(x) 0)
    fidGoodSize[1:TDOnDisk] <- fidOnDisk
    
  } else if (TDOnDisk > TD) {
    warning("Size is bigger than expected, the rest ignored so the size is the same for every fid")
    fidGoodSize <- fidOnDisk(1:TD)
    
  } else {
    fidGoodSize <- fidOnDisk
  }
  
  fidRePart <- fidGoodSize[seq(from = 1, to = TD, by = 2)]
  fidImPart <- fidGoodSize[seq(from = 2, to = TD, by = 2)]
  fid <- complex(real = fidRePart, imaginary = fidImPart)
  
  return(list(fid = fid, params = params))
}
