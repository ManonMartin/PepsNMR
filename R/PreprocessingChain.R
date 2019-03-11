#' @export PreprocessingChain
PreprocessingChain <- function(Fid_data = NULL, Fid_info = NULL, data.path=NULL, readFids = TRUE,
                               groupDelayCorr = TRUE, solventSuppression = TRUE,
                               apodization = TRUE, zerofilling = TRUE,
                               fourierTransform = TRUE, zeroOrderPhaseCorr = TRUE,
                               internalReferencing = TRUE,  baselineCorrection = TRUE,
                               negativeValues0 = TRUE, warping = TRUE, windowSelection = TRUE,
                               bucketing = TRUE,regionRemoval = TRUE, zoneAggregation = TRUE,
                               normalization = TRUE, ...,
                               export = FALSE,
                               format = c("Rdata", "csv","txt"),
                               out.path = '.', filename = "filename", writeArg = c("none", "return", "txt"),
                               verbose = FALSE) {
  
 
  ### Various checks
  
  if (TRUE %in% c(internalReferencing, zeroOrderPhaseCorr, baselineCorrection, negativeValues0,
                  warping, windowSelection, bucketing, regionRemoval, zoneAggregation, normalization) && !fourierTransform ){
    stop("fourierTransform is FALSE but the transformation in the frequency domain is mandatory for some steps.")
  }
  
  
  format <- match.arg(format)
  writeArg <- match.arg(writeArg)
  checkArg(verbose, c("bool"))
  
  # data import checks
  
  if (is.null(Fid_data) + is.null(Fid_info) == 1) {
    stop("Either Fid_data or Fid_info is NULL and the other is not.")
  } else if (is.null(Fid_data) + is.null(Fid_info) == 0) {
    begin_info <- beginTreatment("PreprocessingChain", Fid_data, Fid_info, verbose=verbose)
    Fid_data <- begin_info[["Signal_data"]]
    Fid_info <- begin_info[["Signal_info"]]
    if (readFids) {warning("readFids is TRUE but Fid_data and Fid_info are not NULL !")}
  } else {
    begin_info <- beginTreatment("PreprocessingChain", verbose=verbose)
    if (!readFids) {warning("readFids is FALSE but Fid_data and Fid_info are NULL ! -> no data available")}
    }
  
  
  ### start pre-processing
  
  extraArgs <- list(...) 
  
  if(readFids){
    moreArgs <- list()
    extraArgsNames <- names(extraArgs)[names(extraArgs) %in% methods::formalArgs(ReadFids)]
    moreArgs[extraArgsNames] <- extraArgs[extraArgsNames]
    
    fidList <- do.call(what = "ReadFids", args = c(data.path, moreArgs, verbose=verbose))
    
    data <- fidList[["Fid_data"]]
    Fid_info <- fidList[["Fid_info"]]
  
  } else {
    data <- Fid_data
  }
  

  if(groupDelayCorr){

    moreArgs <- list()
    extraArgsNames <- names(extraArgs)[names(extraArgs) %in% methods::formalArgs(GroupDelayCorrection)]
    moreArgs[extraArgsNames] <- extraArgs[extraArgsNames]
    
    data = do.call(what = "GroupDelayCorrection", args = c(list(Fid_data = data, Fid_info = Fid_info), moreArgs,
                                                            verbose=verbose))

    }
   
  if(solventSuppression){

    moreArgs <- list()
    extraArgsNames <- names(extraArgs)[names(extraArgs) %in% methods::formalArgs(SolventSuppression)]
    moreArgs[extraArgsNames] <- extraArgs[extraArgsNames]
    
    data = do.call(what = "SolventSuppression", args = c(list(Fid_data = data, returnSolvent = FALSE), moreArgs,
                                                         verbose=verbose))

  }

  if(apodization){
    moreArgs <- list()
    extraArgsNames <- names(extraArgs)[names(extraArgs) %in% methods::formalArgs(Apodization)]
    moreArgs[extraArgsNames] <- extraArgs[extraArgsNames]
    
    data = do.call(what = "Apodization", args = c(list(Fid_data = data, Fid_info = Fid_info,
                                                       returnFactor = FALSE), moreArgs,
                                                  verbose=verbose))

  }

  if(zerofilling){
    moreArgs <- list()
    extraArgsNames <- names(extraArgs)[names(extraArgs) %in% methods::formalArgs(ZeroFilling)]
    moreArgs[extraArgsNames] <- extraArgs[extraArgsNames]
    
    data = do.call(what = "ZeroFilling", args = c(list(Fid_data = data), moreArgs, verbose=verbose))
    
  }
  

  if(fourierTransform){
    moreArgs <- list()
    extraArgsNames <- names(extraArgs)[names(extraArgs) %in% methods::formalArgs(FourierTransform)]
    moreArgs[extraArgsNames] <- extraArgs[extraArgsNames]
    
    data = do.call(what = "FourierTransform", args = c(list(Fid_data = data, Fid_info = Fid_info), moreArgs,
                                                       verbose=verbose))

  }


  if(zeroOrderPhaseCorr){
    moreArgs <- list()
    extraArgsNames <- names(extraArgs)[names(extraArgs) %in% methods::formalArgs(ZeroOrderPhaseCorrection)]
    moreArgs[extraArgsNames] <- extraArgs[extraArgsNames]
    
    data = do.call(what = "ZeroOrderPhaseCorrection", args = c(list(Spectrum_data = data,
                                                                    returnAngle = FALSE), moreArgs,
                                                               verbose=verbose))
  }
  
  if(internalReferencing){
    moreArgs <- list()
    extraArgsNames <- names(extraArgs)[names(extraArgs) %in% methods::formalArgs(InternalReferencing)]
    moreArgs[extraArgsNames] <- extraArgs[extraArgsNames]
    data = do.call(what = "InternalReferencing", args = c(list(Spectrum_data = data, 
                                                               Fid_info = Fid_info), moreArgs,
                                                          verbose=verbose))
    
  }
  

  if(baselineCorrection){
    moreArgs <- list()
    extraArgsNames <- names(extraArgs)[names(extraArgs) %in% methods::formalArgs(BaselineCorrection)]
    moreArgs[extraArgsNames] <- extraArgs[extraArgsNames]
    
    data = do.call(what = "BaselineCorrection", args = c(list(Spectrum_data = data,
                                                              returnBaseline = FALSE), moreArgs,
                                                         verbose=verbose))

  }


  if(negativeValues0){
    moreArgs <- list()
    extraArgsNames <- names(extraArgs)[names(extraArgs) %in% methods::formalArgs(NegativeValuesZeroing)]
    moreArgs[extraArgsNames] <- extraArgs[extraArgsNames]
    
    data = do.call(what = "NegativeValuesZeroing", args = c(list(Spectrum_data = data), moreArgs,
                                                            verbose=verbose))

  }

  

  if(warping){
    moreArgs <- list()
    extraArgsNames <- names(extraArgs)[names(extraArgs) %in% methods::formalArgs(Warping)]
    moreArgs[extraArgsNames] <- extraArgs[extraArgsNames]

    data = do.call(what = "Warping", args = c(list(Spectrum_data = data, 
                      returnReference = FALSE, returnWarpFunc = FALSE), moreArgs,
                      verbose=verbose))

  }

  
 
  
  if(windowSelection){
    moreArgs <- list()
    extraArgsNames <- names(extraArgs)[names(extraArgs) %in% methods::formalArgs(WindowSelection)]
    moreArgs[extraArgsNames] <- extraArgs[extraArgsNames]

    data = do.call(what = "WindowSelection", args = c(list(Spectrum_data =  data), moreArgs,
                                                      verbose=verbose))

  }


  if(bucketing){
    moreArgs <- list()
    extraArgsNames <- names(extraArgs)[names(extraArgs) %in% methods::formalArgs(Bucketing)]
    moreArgs[extraArgsNames] <- extraArgs[extraArgsNames]

    data = do.call(what = "Bucketing", args = c(list(Spectrum_data = data), moreArgs,
                                                verbose=verbose))

  }


  if(regionRemoval){
    moreArgs <- list()
    extraArgsNames <- names(extraArgs)[names(extraArgs) %in% methods::formalArgs(RegionRemoval)]
    moreArgs[extraArgsNames] <- extraArgs[extraArgsNames]

    data = do.call(what = "RegionRemoval", args = c(list(Spectrum_data = data), moreArgs,
                                                    verbose=verbose))

  }



  if(zoneAggregation){
    moreArgs <- list()
    extraArgsNames <- names(extraArgs)[names(extraArgs) %in% methods::formalArgs(ZoneAggregation)]
    moreArgs[extraArgsNames] <- extraArgs[extraArgsNames]

    data = do.call(what = "ZoneAggregation", args = c(list(Spectrum_data = data), moreArgs,
                                                      verbose=verbose))

  }



  if(normalization){
    moreArgs <- list()
    extraArgsNames <- names(extraArgs)[names(extraArgs) %in% methods::formalArgs(Normalization)]
    moreArgs[extraArgsNames] <- extraArgs[extraArgsNames]

    data = do.call(what = "Normalization", args = c(list(Spectrum_data = data, returnFactor = FALSE), moreArgs,
                                                    verbose=verbose))

  }

  Spectra = data
  
  if (export) {

    if (format == "Rdata"){
      save(Spectra, Fid_info, file = file.path(out.path, paste0(filename, ".RData")))
    } else if (format == "csv"){
      utils::write.csv(Spectra, file = file.path(out.path, paste0(filename, "_SpectralIntensities.csv")))
      utils::write.csv(Fid_info, file = file.path(out.path, paste0(filename, "_FidInfo.csv")))
    } else { # txt
      utils::write.table(Spectra, file = file.path(out.path, paste0(filename, "_SpectralIntensities.txt")),
                  sep="\t",row.names=TRUE)
      utils::write.table(Fid_info, file = file.path(out.path, paste0(filename, "_FidInfo.txt")),
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
  return(list(Spectrum_data = endTreatment("PreprocessingChain", begin_info, data, verbose=verbose),
              Fid_info = Fid_info, arguments = arguments))
  } else {return(list(Spectrum_data = endTreatment("PreprocessingChain", begin_info, data, verbose=verbose),
                    Fid_info = Fid_info))}
  
  
}