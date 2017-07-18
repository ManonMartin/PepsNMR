#' @export PreprocessingChain
PreprocessingChain <- function(Fid_data = NULL, Fid_info = NULL, data.path=NULL, readFids = TRUE,
                               groupDelayCorr = TRUE, solventSuppression = TRUE,
                               apodization = TRUE,
                               fourierTransform = TRUE, internalReferencing = TRUE,
                               zeroOrderPhaseCorr = TRUE,  baselineCorrection = TRUE,
                               negativeValues0 = TRUE, warping = TRUE, windowSelection = TRUE,
                               bucketing = TRUE,regionRemoval = TRUE, zoneAggregation = TRUE,
                               normalization = TRUE, ...,
                               export = FALSE,
                               format = c("Rdata", "csv","txt"),
                               out.path = '.', filename = "filename", writeArg = c("none", "return", "txt")) {
  
 
  ### Various checks
  
  if (TRUE %in% c(internalReferencing, zeroOrderPhaseCorr, baselineCorrection, negativeValues0,
                  warping, windowSelection, bucketing, regionRemoval, zoneAggregation, normalization) && !fourierTransform ){
    stop("fourierTransform is FALSE but the transformation in the frequency domain is mandatory for some steps.")
  }
  
  
  format <- match.arg(format)
  writeArg <- match.arg(writeArg)
  
  
  # data import checks
  
  if (is.null(Fid_data) + is.null(Fid_info) == 1) {
    error("Either Fid_data or Fid_info is NULL and the other is not.")
  } else if (is.null(Fid_data) + is.null(Fid_info) == 0) {
    begin_info <- beginTreatment("PreprocessingChain", Fid_data, Fid_info)
    Fid_data <- begin_info[["Signal_data"]]
    Fid_info <- begin_info[["Signal_info"]]
    if (readFids) {warning("readFids is TRUE but Fid_data and Fid_info are not NULL !")}
  } else {
    begin_info <- beginTreatment("PreprocessingChain")
    if (!readFids) {warning("readFids is FALSE but Fid_data and Fid_info are NULL ! -> no data available")}
    }
  
  
  ### start pre-processing
  
  extraArgs <- list(...) 
  
  if(readFids){
    moreArgs <- list()
    extraArgsNames <- names(extraArgs)[names(extraArgs) %in% formalArgs(ReadFids)]
    moreArgs[extraArgsNames] <- extraArgs[extraArgsNames]
    
    fidList <- do.call(what = "ReadFids", args = c(data.path, moreArgs))
    
    data <- fidList[["Fid_data"]]
    Fid_info <- fidList[["Fid_info"]]
  
  } else {
    data <- Fid_data
  }
  

  if(groupDelayCorr){

    moreArgs <- list()
    extraArgsNames <- names(extraArgs)[names(extraArgs) %in% formalArgs(GroupDelayCorrection)]
    moreArgs[extraArgsNames] <- extraArgs[extraArgsNames]
    
    data = do.call(what = "GroupDelayCorrection", args = c(list(Fid_data = data, Fid_info = Fid_info), moreArgs))

    }
   
  if(solventSuppression){

    moreArgs <- list()
    extraArgsNames <- names(extraArgs)[names(extraArgs) %in% formalArgs(SolventSuppression)]
    moreArgs[extraArgsNames] <- extraArgs[extraArgsNames]
    
    data = do.call(what = "SolventSuppression", args = c(list(Fid_data = data, returnSolvent = F), moreArgs))

  }

  if(apodization){
    moreArgs <- list()
    extraArgsNames <- names(extraArgs)[names(extraArgs) %in% formalArgs(Apodization)]
    moreArgs[extraArgsNames] <- extraArgs[extraArgsNames]
    
    data = do.call(what = "Apodization", args = c(list(Fid_data = data, Fid_info = Fid_info,
                                                       returnFactor = F), moreArgs))

  }


  if(fourierTransform){
    moreArgs <- list()
    extraArgsNames <- names(extraArgs)[names(extraArgs) %in% formalArgs(FourierTransform)]
    moreArgs[extraArgsNames] <- extraArgs[extraArgsNames]
    
    data = do.call(what = "FourierTransform", args = c(list(Fid_data = data, Fid_info = Fid_info), moreArgs))

  }


  if(internalReferencing){
    moreArgs <- list()
    extraArgsNames <- names(extraArgs)[names(extraArgs) %in% formalArgs(InternalReferencing)]
    moreArgs[extraArgsNames] <- extraArgs[extraArgsNames]
    data = do.call(what = "InternalReferencing", args = c(list(Spectrum_data = data, 
                                                               Fid_info = Fid_info), moreArgs))

  }


  if(zeroOrderPhaseCorr){
    moreArgs <- list()
    extraArgsNames <- names(extraArgs)[names(extraArgs) %in% formalArgs(ZeroOrderPhaseCorrection)]
    moreArgs[extraArgsNames] <- extraArgs[extraArgsNames]
    
    data = do.call(what = "ZeroOrderPhaseCorrection", args = c(list(Spectrum_data = data,
                                                                    returnAngle = FALSE), moreArgs))
  }


  if(baselineCorrection){
    moreArgs <- list()
    extraArgsNames <- names(extraArgs)[names(extraArgs) %in% formalArgs(BaselineCorrection)]
    moreArgs[extraArgsNames] <- extraArgs[extraArgsNames]
    
    data = do.call(what = "BaselineCorrection", args = c(list(Spectrum_data = data,
                                                              returnBaseline = F), moreArgs))

  }


  if(negativeValues0){
    moreArgs <- list()
    extraArgsNames <- names(extraArgs)[names(extraArgs) %in% formalArgs(NegativeValuesZeroing)]
    moreArgs[extraArgsNames] <- extraArgs[extraArgsNames]
    
    data = do.call(what = "NegativeValuesZeroing", args = c(list(Spectrum_data = data), moreArgs))

  }

  

  if(warping){
    moreArgs <- list()
    extraArgsNames <- names(extraArgs)[names(extraArgs) %in% formalArgs(Warping)]
    moreArgs[extraArgsNames] <- extraArgs[extraArgsNames]

    data = do.call(what = "Warping", args = c(list(Spectrum_data = data, 
                      returnReference = FALSE, returnWarpFunc = FALSE), moreArgs))

  }

  
 
  
  if(windowSelection){
    moreArgs <- list()
    extraArgsNames <- names(extraArgs)[names(extraArgs) %in% formalArgs(WindowSelection)]
    moreArgs[extraArgsNames] <- extraArgs[extraArgsNames]

    data = do.call(what = "WindowSelection", args = c(list(Spectrum_data =  data), moreArgs))

  }


  if(bucketing){
    moreArgs <- list()
    extraArgsNames <- names(extraArgs)[names(extraArgs) %in% formalArgs(Bucketing)]
    moreArgs[extraArgsNames] <- extraArgs[extraArgsNames]

    data = do.call(what = "Bucketing", args = c(list(Spectrum_data = data), moreArgs))

  }


  if(regionRemoval){
    moreArgs <- list()
    extraArgsNames <- names(extraArgs)[names(extraArgs) %in% formalArgs(RegionRemoval)]
    moreArgs[extraArgsNames] <- extraArgs[extraArgsNames]

    data = do.call(what = "RegionRemoval", args = c(list(Spectrum_data = data), moreArgs))

  }



  if(zoneAggregation){
    moreArgs <- list()
    extraArgsNames <- names(extraArgs)[names(extraArgs) %in% formalArgs(ZoneAggregation)]
    moreArgs[extraArgsNames] <- extraArgs[extraArgsNames]

    data = do.call(what = "ZoneAggregation", args = c(list(Spectrum_data = data), moreArgs))

  }



  if(normalization){
    moreArgs <- list()
    extraArgsNames <- names(extraArgs)[names(extraArgs) %in% formalArgs(Normalization)]
    moreArgs[extraArgsNames] <- extraArgs[extraArgsNames]

    data = do.call(what = "Normalization", args = c(list(Spectrum_data = data, returnFactor = F), moreArgs))

  }

  Spectra = data
  
  if (export) {

    if (format == "Rdata"){
      save(Spectra, Fid_info, file = file.path(out.path, paste0(filename, ".RData")))
    } else if (format == "csv"){
      write.csv(Spectra, file = file.path(out.path, paste0(filename, "_SpectralIntensities.csv")))
      write.csv(Fid_info, file = file.path(out.path, paste0(filename, "_FidInfo.csv")))
    } else { # txt
      write.table(Spectra, file = file.path(out.path, paste0(filename, "_SpectralIntensities.txt")),
                  sep="\t",row.names=TRUE)
      write.table(Fid_info, file = file.path(out.path, paste0(filename, "_FidInfo.txt")),
                  sep="\t",row.names=TRUE)
    }

  }
  
  # write lines in txt or return a list with all the arguments
  arguments <- c(as.list(environment()), list(...))
  arguments[names(arguments) %in% c("Fid_info", "data", "fidList", "Spectra",
                                     "extraArgsNames", "moreArgs","extraArgs",
                                     "begin_info","Fid_data" )] = NULL

  if (writeArg=="txt") {
    sink(file.path(out.path,"PreprocessingChainArguments.txt"))
    print(arguments)
    sink()
  }

  
  if (writeArg=="return") {
  return(list(Spectrum_data = endTreatment("PreprocessingChain", begin_info, data), Fid_info = Fid_info, arguments = arguments))
} else {return(list(Spectrum_data = endTreatment("PreprocessingChain", begin_info, data), Fid_info = Fid_info))}
  
  
}