##################################################
##################################################

###  H-NMR spectra pre-treatment with SOAP    ###

##################################################
##################################################


###################
# INITIALISATION
###################

## Arguments
###################


# dataname : includes the dataset name and additional infos
# save : if the data needs to be saved in a Rdata file 
# saveall : if the data from each step need to be saved in a Rdata file 
# data.path : path to the main directory
# out.path : output path for datasets and graphs
# nspectr : choice of the observation chosen for the graphs
#  If TRUE, 

## TRUE or FALSE :
# Fopc : First order phase correction
# Ss : Solvent suppression
# A : Apodization
# Zopc = Zero order phase correction
# Bc : Baseline correction
# Zsnv : Zero setting of negative values
# W : Warping
# B : Bucketing
# Zs : Zone suppression
# Za :  Zone aggregation
# N : Normalisation
# ImpG : impression graphique if TRUE
# RetArgs : save function Arguments 
  
## Default values parameters for the inner functions
# SolventSuppression : lambda.Ss = 1e6
# Apodization : type.A = "exp"
# BaselineCorrection : p.Bc = 0.05, lambda.Bc = 1e7
# Warping : normalization.type.W ="median" , reference.choosing.W="fixed" , reference.W= 1 , K.W=3, L.W=40,lambda.smooth.W=0, deg.W=3, lambda.bspline.W=0.01, kappa.W=0.0001
# PPMConversion  : shiftHandling.Pc ="cut"
# WindowSelection  : from.Ws=0.2, to.Ws =10
# Bucketing : m.B=500
# RegionRemoval : 
#   -for serum = fromto.Zs=list(Water =c(4.5, 5.04), Lactate=c(1.36, 2.8))
#   -for urine = fromto.Zs=list(Water =c(4.5, 5), Uree=c(4.5, 6))
# ZoneAggregation : fromto.Za=list(Citrate =c(2.5, 2.7))
# Normalization : type.N = "mean"


## VALUES
##############
# Spectra : returns the pretreated spectra
# List spectra : list of spectra after each step
# graphs


##################################################
##################################################


#' @export PreprocessingChain

PreprocessingChain = function(dataname = "Dataset", data.path = getwd(), out.path = getwd(), 
                        nspectr = 1, save = FALSE, saveall = FALSE, ImpG= FALSE, RetArgs = TRUE, export = NULL, 
                        Fopc = TRUE, Ss = TRUE, A = TRUE, Zopc = TRUE, Bc = TRUE, 
                        Zsnv = TRUE, W = TRUE, B = TRUE, Zs = TRUE, Za=FALSE, N = TRUE,
                        l=1, subdirs = FALSE, #ReadFids
                        group_delay=NULL, # FOPC
                        lambda.ss=1e6, ptw.ss=TRUE,  # SolventSuppression
                        DT=NULL,type.apod = "exp",phase=0, rectRatio=1/2, gaussLB=1, expLB=1, # Apodization
                        SW_h=NULL, # FourierTransform
                        Angle = NULL,  p.zo=0.8, # ZeroOrderPhaseCorrection
                        ptw.bc=TRUE, maxIter = 42,lambda.bc=1e7, p.bc=0.05, eps=1e-8, # BaselineCorrection
                        shiftHandling="cut", c = 2, #PPMConversion
                        
                        normalization.type="median", from.normW=3.05, to.normW=4.05,reference.choosing="fixed", 
                        reference=1,optim.crit="RMS", ptw.wp=F, K=3, L=40,
                        lambda.smooth=0, deg=3, lambda.bspline=0.01, kappa=0.0001,
                        max_it_Bspline=10, returnReference=F,  #Warping
                        from.ws = 0.2, to.ws = 10, reverse.axis = TRUE, # WindowSelection
                        m = 500, # Bucketing
                        typeofspectra = "serum",type.rr =  "zero", fromto.rr=list(Water =c(4.5, 5.1)), # RegionRemoval
                        fromto.za = list(Citrate =c(2.5, 2.7)), # ZoneAggregation
                        type.norm="mean", from.norm=3.05, to.norm=4.05, ref.norm=1 # Normalization
                        ) { 


 
#_______________________________________________________
#############################
#############################
#   PRETREATMENT WORKFLOW
#############################
#############################
#_______________________________________________________

  
#_______________________________________________________
  
cat("##############################################","\n")
cat("PRETREATMENT OF ", dataname,"\n")
cat("############################################## \n")
  
if (saveall == TRUE) {
  ll = length(which(c(Fopc , Ss , A , Zopc, Bc ,Zsnv , W , B, Zs , Za, N)==TRUE)) # optional steps
  len = ll + 4 # nombre d'etapes totales
    
  PretreatedSpectrabyStep = vector("list", len) # creation d une liste pour sauver les donnees de chaque etape
  step=1
}
  #_______________________________________________________  



##########################
## Load FIDs and info
##########################

fidList <-ReadFids(data.path, l=l, subdirs = subdirs)

# Initial dataset 
Fid_info <- fidList[["Fid_info"]]
Fid_data <- fidList[["Fid_data"]]
Fid_data0 = Fid_dataB = Fid_data


if (saveall == TRUE) {
  PretreatedSpectrabyStep[[step]] = Fid_data0 
  names(PretreatedSpectrabyStep)[step] <- "InitialFID.Data"
  step = step + 1
}
  

  
################################
## FirstOrderPhaseCorrection
################################

if (Fopc ==  TRUE ){

  Fid_data <- FirstOrderPhaseCorrection(Fid_data, Fid_info = Fid_info, group_delay=group_delay)
  
  if (saveall == TRUE) {
    PretreatedSpectrabyStep[[step]] = Fid_data 
    names(PretreatedSpectrabyStep)[step] <- "FirstOrderPhaseCorrection.Data"
    step = step + 1
  }
 

  if (ImpG==TRUE) {
    grDevices::pdf(paste0(out.path,"/FirstOrderPhCorrectREAL.pdf"),width=13, height=8)
    graphics::par(mfrow=c(2,2))
    graphics::plot(Re(Fid_dataB[nspectr,]), col="blue", type="l", ylab = "Intensity", xlab="Time", main="Original FID \n Real part")
    graphics::plot(Re(Fid_dataB[nspectr,1:200]), col="blue", ylab = "Intensity",xlab="Time", type="l", main="Original FID (zoom)\n Real part ")
    graphics::plot(Re(Fid_data[nspectr,]), col="blue",  ylab = "Intensity",xlab="Time", type="l", main="1st Order Phase Correction \n Real part")
    graphics::plot(Re(Fid_data[nspectr,1:200]), col="blue", ylab = "Intensity",xlab="Time", type="l", main="1st Order Phase Correction (zoom)\n Real part")
    grDevices::dev.off()
    
    grDevices::pdf(paste0(out.path,"/FirstOrderPhCorrectIMAG.pdf"),width=13, height=8)
    graphics::par(mfrow=c(2,2))
    graphics::plot(Im(Fid_dataB[nspectr,]), col="blue", ylab = "Intensity",xlab="Time", type="l", main="Original FID \n Imaginary part")
    graphics::plot(Im(Fid_dataB[nspectr,1:200]), col="blue", ylab = "Intensity",xlab="Time", type="l", main="Original FID (zoom)\n Imaginary part ")
    graphics::plot(Im(Fid_data[nspectr,]), col="blue", ylab = "Intensity",xlab="Time", type="l", main="1st Order Phase Correction \n Imaginary part")
    graphics::plot(Im(Fid_data[nspectr,1:200]), col="blue", ylab = "Intensity",xlab="Time", type="l", main="1st Order Phase Correction (zoom)\n Imaginary part ")
    grDevices::dev.off()
  }
}




##########################
# SolventSuppression
##########################

if (Ss ==  TRUE ){

  Fid_dataB = Fid_data 
  Ss.res = SolventSuppression(Fid_data,returnSolvent=TRUE, lambda.ss=lambda.ss, ptw.ss=ptw.ss, plotSolvent=FALSE)
  Fid_data = Ss.res[["Fid_data"]]
  SolventRe = Ss.res[["SolventRe"]]
  
  if (saveall == TRUE) {
    PretreatedSpectrabyStep[[step]] = Fid_data
    names(PretreatedSpectrabyStep)[step] <- "SolventSuppression.Data"
    step = step + 1
  }

  if (ImpG==TRUE) {
    grDevices::pdf(paste0(out.path,"/SolventSuppression1.pdf"),width=8, height=4)
    graphics::plot(Re(Fid_data1[nspectr,]), ylim = c(min(Re(Fid_data1[nspectr,])), max(Re(Fid_data1[nspectr,]))), col="blue", type="l", ylab = "Intensity", xlab="Time", main="Real part of the FID and solvent residuals signal")
    graphics::lines(SolventRe[nspectr,],col="red" )
    graphics::legend("topright", legend = "Solvent signal", col="red",  lty = 1)
    grDevices::dev.off()
    
    grDevices::pdf(paste0(out.path,"/SolventSuppression2.pdf"),width=13, height=8)
    graphics::par(mfrow=c(2,2))
    graphics::plot(Re(Fid_dataB[nspectr,1:100]), col="blue", type="l", ylab = "Intensity", xlab="Time", main="Before solvent suppression (zoom)\n Real part ")
    graphics::plot(Im(Fid_dataB[nspectr,1:100]), col="blue", type="l",  ylab = "Intensity", xlab="Time", main="Before solvent suppression (zoom)\n Imaginary part ")
    graphics::plot(Re(Fid_data[nspectr,1:100]), col="blue", type="l",  ylab = "Intensity", xlab="Time", main="After solvent suppression (zoom)\n Real part ")
    graphics::plot(Im(Fid_data[nspectr,1:100]), col="blue", type="l",  ylab = "Intensity", xlab="Time", main="After solvent suppression (zoom)\n Imaginary part ")
    grDevices::dev.off()
  }

}



##########################
# Apodization
##########################
if (A ==  TRUE ){
  
  Fid_dataB = Fid_data 
  Apod.res = Apodization(Fid_data, Fid_info = Fid_info, returnFactor=T, DT=DT,
                         type.apod = type.apod,phase=phase, rectRatio=rectRatio, gaussLB=gaussLB, 
                         expLB=expLB, plotWindow=FALSE)
  Fid_data = Apod.res[["Fid_data"]]
  ApodFactor = Apod.res[["Factor"]]
  


  if (saveall == TRUE) {
    PretreatedSpectrabyStep[[step]] = Fid_data 
    names(PretreatedSpectrabyStep)[step] <- "Apodization.Data"
    step = step + 1
  }
  
  if (ImpG==TRUE) {
    grDevices::pdf(paste0(out.path,"/Apodization1.pdf"),width=13, height=8)
    graphics::par(mfrow=c(2,2))
    graphics::plot(Re(Fid_dataB[nspectr,]), col="blue", type="l", xlab="Time", main="Before apodization \n Real part")
    graphics::plot(Im(Fid_dataB[nspectr,]), col="blue", type="l", xlab="Time", main="Before apodization \n Imaginary part")
    graphics::plot(Re(Fid_data[nspectr,]), col="blue", type="l", xlab="Time", main="After apodization \n Real part")
    graphics::plot(Im(Fid_data[nspectr,]), col="blue", type="l", xlab="Time", main="After apodization \n Imaginary part")
    grDevices::dev.off()
    
    grDevices::pdf(paste0(out.path,"/Apodization2.pdf"),width=8, height=4)
    graphics::plot(Apod.res[["Factor"]], col="red",type="l", main = "Apodization factor")
    grDevices::dev.off()
  }

}

##########################
# FourierTransform
##########################


RawSpect_data = FourierTransform(Fid_data, Fid_info = Fid_info, SW_h = SW_h)

FourierTransformData = Re(RawSpect_data)

if (save == TRUE) {
  save(FourierTransformData = FourierTransformData,  file=paste0(out.path, "/",dataname ,"FourierTransformData.RData"))
}

if (saveall == TRUE) {
  PretreatedSpectrabyStep[[step]] = RawSpect_data
  names(PretreatedSpectrabyStep)[step] <- "FourierTransform.Data"
  step = step + 1
}


if (ImpG==TRUE) {
  grDevices::pdf(paste0(out.path,"/FourierTransform.pdf"),width=13, height=14)
  graphics::par(mfrow=c(3,1))
  aa = 0.2*length(RawSpect_data[nspectr,])
  bb = 0.6*length(RawSpect_data[nspectr,])
  cc = 0.06*length(RawSpect_data[nspectr,])
  
  graphics::plot(Re(RawSpect_data[nspectr,]), col="blue", type="l", xlab="Frequency", main="After Fourier Transform \n Real part")
  graphics::plot(Re(RawSpect_data[nspectr,aa:bb]), col="blue", xaxt = "n", type="l", xlab="Frequency", main="After Fourier Transform \n Real part (zoom)")
  at=seq(1,(bb-aa), cc)
  graphics::axis(side=1, at=at, labels  = c(aa:bb)[at])
  graphics::plot(Im(RawSpect_data[nspectr,aa:bb]), col="blue", xaxt = "n", type="l", xlab="Frequency", main="After Fourier Transform \n Imaginary part (zoom)")
  at=seq(1,(bb-aa), cc)
  graphics::axis(side=1, at=at, labels = c(aa:bb)[at])
  grDevices::dev.off()
  ###
}

##########################
# ZeroOrderPhaseCorrection
##########################


if (Zopc ==  TRUE ){

  RawSpect_dataB = RawSpect_data
  RawSpect_data = ZeroOrderPhaseCorrection(RawSpect_data, plot_rms = NULL, returnAngle = FALSE, createWindow = FALSE, Angle = Angle,   p.zo=p.zo, plot_spectra = FALSE)

  if (saveall == TRUE) {
    PretreatedSpectrabyStep[[step]] = RawSpect_data
    names(PretreatedSpectrabyStep)[step] <- "ZeroOrderPhaseCorrection.Data"
    step = step + 1
  }
  
  if (ImpG==TRUE) {
    grDevices::pdf(paste0(out.path,"/ZeroOrderPhaseCorrection.pdf"),width=13, height=14)
    graphics::par(mfrow=c(3,1))
    
    aa = 0.2*length(RawSpect_dataB[nspectr,])
    bb = 0.6*length(RawSpect_dataB[nspectr,])
    cc = 0.06*length(RawSpect_dataB[nspectr,])
    
    graphics::plot(Re(RawSpect_dataB[nspectr,aa:bb]), col="blue",  xaxt = "n", type="l", xlab="Frequency", main="Before 0 Order Phase Correction \n Real part (zoom)")
    at=seq(1,(bb-aa), cc)
    graphics::axis(side=1, at=at, labels = c(aa:bb)[at])
    graphics::plot(Re(RawSpect_data[nspectr,aa:bb]), col="blue",  xaxt = "n", type="l", xlab="Frequency", main="After 0 Order Phase Correction \n Real part (zoom)")
    at=seq(1,(bb-aa), cc)
    graphics::axis(side=1, at=at, labels = c(aa:bb)[at])
    graphics::plot(Im(RawSpect_data[nspectr,aa:bb]), col="blue",  xaxt = "n", type="l", xlab="Frequency", main="After 0 Order Phase Correction \n Imaginary part (zoom)")
    at=seq(1,(bb-aa), cc)
    graphics::axis(side=1, at=at, labels = c(aa:bb)[at])
    grDevices::dev.off()
    }

}

##########################
# BaselineCorrection
##########################
if (Bc ==  TRUE ){

  RawSpect_dataB = RawSpect_data
  BC.res =  BaselineCorrection(RawSpect_data, returnBaseline=TRUE, lambda.bc=lambda.bc, p.bc=p.bc, eps=eps)
  baseline = BC.res[["Baseline"]]
  RawSpect_data = BC.res[["RawSpect_data"]]



  if (saveall == TRUE) {
    PretreatedSpectrabyStep[[step]] = RawSpect_data
    names(PretreatedSpectrabyStep)[step] <- "BaselineCorrection.Data"
    step = step + 1
  }

  if (ImpG==TRUE) {
    grDevices::pdf(paste0(out.path,"/BaselineCorrection.pdf"),width=13, height=8)
    graphics::par(mfrow=c(2,2))
    
    aa = 0.42*length(RawSpect_dataB[nspectr,])
    bb = 0.46*length(RawSpect_dataB[nspectr,])
    
    graphics::plot(Re(RawSpect_dataB[nspectr,]), col="blue", type="l", ylab = "Intensity", xlab="Frequency", main="Before Baseline Correction")
    graphics::abline(h=0, lty = 2)
    graphics::lines(baseline[nspectr,], col="red")
    
    graphics::plot(Re(RawSpect_dataB[nspectr,aa:bb]), col="blue", xaxt="n", type="l", ylab = "Intensity", xlab="Frequency", main="Before Baseline Correction  \n zoom")
    graphics::axis(side = 1, at = seq(0,(bb-aa), 100), labels = seq(aa,bb, 100))
    graphics::abline(h=0, lty = 2)
    graphics::lines(baseline[nspectr,aa:bb], col="red")
    
    graphics::plot(Re(RawSpect_data[nspectr,]), col="blue", type="l", ylab = "Intensity", xlab="Frequency", main="After BaselineCorrection")
    graphics::abline(h=0, lty = 2)
    
    graphics::plot(Re(RawSpect_data[nspectr,aa:bb]),  col="blue", xaxt="n", type="l", ylab = "Intensity", xlab="Frequency", main="After BaselineCorrection \n zoom")
    graphics::axis(side = 1, at = seq(0,(bb-aa), 100), labels = seq(aa,bb, 100))
    graphics::abline(h=0, lty = 2)
    
    grDevices::dev.off()
  }

}

##########################
# NegativeValuesZeroing
##########################

if (Zsnv ==  TRUE ){
  
  RawSpect_dataB = RawSpect_data
  RawSpect_data = NegativeValuesZeroing(RawSpect_data)

  if (saveall == TRUE) {
    PretreatedSpectrabyStep[[step]] = RawSpect_data
    names(PretreatedSpectrabyStep)[step] <- "NegativeValuesZeroing.Data"
    step = step + 1
  }
  
  if (ImpG==TRUE) {
    grDevices::pdf(paste0(out.path,"/NegativeValuesZeroing.pdf"), width = 13, height = 8)
    graphics::par(mfrow=c(2,1))
    
    aa = 0.2*length(RawSpect_dataB[nspectr,])
    bb = 0.6*length(RawSpect_dataB[nspectr,])
    cc = 0.06*length(RawSpect_dataB[nspectr,])
    
    graphics::plot(Re(RawSpect_dataB[nspectr,aa:bb]), col="blue", xaxt = "n",type="l", xlab="Frequency", main="Before Negative Values Zeroing \n Real part (zoom)")
    at=seq(1,(bb-aa), cc)
    graphics::axis(side=1, at=at, labels = c(aa:bb)[at])
    graphics::plot(Re(RawSpect_data[nspectr,aa:b]), col="blue", xaxt = "n",type="l", xlab="Frequency", main="After Negative Values Zeroing \n Real part (zoom)")
    at=seq(1,(bb-aa), cc)
    graphics::axis(side=1, at=at, labels = c(aa:bb)[at])
    grDevices::dev.off()
  }
  
}


##########################
# PPMConversion 
##########################

Spectrum_data = PPMConversion(RawSpect_data, RawSpect_info = Fid_info, shiftHandling=shiftHandling, c = c)


if (saveall == TRUE) {
  PretreatedSpectrabyStep[[step]] = Spectrum_data
  names(PretreatedSpectrabyStep)[step] <- "PPMConversion.Data"
  step = step + 1
}

if (ImpG==TRUE) {
  grDevices::pdf(paste0(out.path,"/PPMConversion.pdf"),width=13, height=8)
  graphics::par(mfrow=c(2,1))
  graphics::plot(Re(RawSpect_data[nspectr,]), col="blue", type="l", xlab="Frequency", main="Before PPM Conversion \n Real part")
  xat=round(as.numeric(colnames(Spectrum_data)),5)
  graphics::plot(xat,Re(Spectrum_data[nspectr,]), col="blue", type="l", xlab="ppm", main="After PPM Conversion \n Real part")
  grDevices::dev.off()
}


############
# Warping
############

if (W ==  TRUE ){
  Spectrum_dataB = Spectrum_data
  Spectrum_data  = Warping(RawSpect_data=Spectrum_data, normalization.type=normalization.type, 
                      from.normW=from.normW, to.normW=to.normW,reference.choosing=reference.choosing, 
                      reference=reference,optim.crit=optim.crit, ptw.wp=ptw.wp, K=K, L=L,
                      lambda.smooth=lambda.smooth, deg=deg, lambda.bspline=lambda.bspline, kappa=kappa,
                      max_it_Bspline=max_it_Bspline, returnReference=returnReference)
  
  
  if (saveall == TRUE) {
    PretreatedSpectrabyStep[[step]] = Spectrum_data
    names(PretreatedSpectrabyStep)[step] <- "Warping.Data"
    step = step + 1
  }
  
  if (ImpG==TRUE) {
    
    grDevices::pdf(paste0(out.path,"/Warping.pdf"), width = 13, height = 10)
    graphics::par(mfrow=c(2,1),mar=c(4.1, 4.1, 4.1, 9.5), xpd=TRUE)
    
    if (dim(Spectrum_dataB)[1]>=3) {
      f=sample(c(1:dim(Spectrum_dataB)[1]), 3)
    }else {f = 1}

    
    aa = 0.4*length(RawSpect_dataB[f[1],])
    bb = 0.47*length(RawSpect_dataB[f[1],])
    cc = 0.015*length(RawSpect_dataB[f[1],])
    
    graphics::plot(Re(Spectrum_dataB[nspectr,aa:bb]),  col="red", xaxt = "n",ylab = "Intensity",ylim=c(0, max(Re(Spectrum_dataB[c(nspectr, f),aa:bb]))), type="l", xlab="Frequency", main="Before Warping (zoom)")
    graphics::axis(side = 1, at = seq(0,(bb-aa), cc), labels = seq(aa,bb, cc))   
    graphics::legend("topright", inset=c(-0.22,-0.12), legend=c("Ref. spectrum", "Warped  spectra"), lty = 1, cex=0.8, col=c("red", "blue"))    
    for (j in f) {
      graphics::lines(Re(Spectrum_dataB[j,aa:b]), col="blue", type="l")}
    graphics::grid(20, NA, lty = 1, lwd = 0.7, col="gray44")
    
    graphics::plot(Re(Spectrum_data[nspectr,aa:b]), col="red", xaxt = "n",ylab = "Intensity",ylim=c(0, max(Re(Spectrum_data[c(nspectr, f),aa:bb]))), type="l", xlab="Frequency", main="After Warping (zoom)")
    graphics::axis(side = 1, at = seq(0,(bb-aa), cc), labels = seq(aa,bb, cc))    
    graphics::legend("topright", inset=c(-0.22,-0.12), legend=c("Ref. spectrum", "Warped  spectra"), lty = 1, cex=0.8, col=c("red", "blue"))    
    for (j in f) {
      graphics::lines(Re(Spectrum_data[j,aa:bb]), type="l", col="blue")}
    graphics::grid(20, NA, lty = 1, lwd = 0.7, col="gray44")
    grDevices::dev.off()
  }
  
}


##########################
# WindowSelection
##########################
Spectrum_dataB = Spectrum_data
Spectrum_data = WindowSelection(Spectrum_data, from.ws = from.ws, to.ws = to.ws, reverse.axis = reverse.axis)


if (saveall == TRUE) {
  PretreatedSpectrabyStep[[step]] = Spectrum_data
  names(PretreatedSpectrabyStep)[step] <- "WindowSelection.Data"
  step = step + 1
}


if (ImpG==TRUE) {
  grDevices::pdf(paste0(out.path,"/WindowSelection.pdf"), width = 10, height = 9)
  graphics::par(mfrow=c(2,1), mar=c(4,4,2,2))
  cc = 0.03*length(RawSpect_dataB[nspectr,])
  xat=round(as.numeric(colnames(Spectrum_dataB)),5)
  graphics::plot(xat,Re(Spectrum_dataB[nspectr,]),col="blue",  type="l",  ylab = "Intensity", xlab="ppm", main="Before Window Selection")
  xat=round(as.numeric(colnames(Spectrum_data)),5)
  graphics::plot( Re(Spectrum_data[nspectr,]), col="blue",  ylab = "Intensity", xaxt="n", type="l", xlab="ppm", main="After Window Selection")
  at=seq(1,length(xat), cc)
  graphics::axis(side=1, at=at, labels=round(xat[at],2))
  grDevices::dev.off()
}
 
################ 
# Bucketing
################

if (B ==  TRUE ){
  Spectrum_dataB = Spectrum_data
  Spectrum_data = Bucketing(Spectrum_data, m = m)

  if (saveall == TRUE) {
    PretreatedSpectrabyStep[[step]] = Spectrum_data
    names(PretreatedSpectrabyStep)[step] <- "Bucketing.Data"
    step = step + 1
  }
  
  if (ImpG==TRUE) {
    grDevices::pdf(paste0(out.path,"/Bucketing.pdf"), width = 10, height = 9)
    graphics::par(mfrow=c(2,1), mar=c(4,4,2,2))
    cc = 0.03*length(RawSpect_dataB[nspectr,])
    xat=round(as.numeric(colnames(Spectrum_dataB)),5)
    graphics::plot(Re(Spectrum_dataB[nspectr,]), col="blue",ylab = "Intensity", xaxt="n", type="l", xlab="ppm", main="Before Bucketing")
    at=seq(1,length(xat), cc)
    graphics::axis(side=1, at=at, labels=round(xat[at],2))
    xat=round(as.numeric(colnames(Spectrum_data)),5)
    graphics::plot(Re(Spectrum_data[nspectr,]), col="blue",ylab = "Intensity", xaxt="n", type="l", xlab="ppm", main="After Bucketing")
    at=seq(1,length(xat), 30)
    graphics::axis(side=1, at=at, labels=round(xat[at],2))
    grDevices::dev.off()
  }
  
  
}

##########################
# RegionRemoval
##########################

if (Zs ==  TRUE ){
  Spectrum_dataB = Spectrum_data
  Spectrum_data = RegionRemoval(Spectrum_data,typeofspectra = typeofspectra, 
                                type.rr = type.rr, fromto.rr=fromto.rr) 


  if (saveall == TRUE) {
    PretreatedSpectrabyStep[[step]] = Spectrum_data
    names(PretreatedSpectrabyStep)[step] <- "RegionRemoval.Data"
    step = step + 1
  }
  
  if (ImpG==TRUE) {
    grDevices::pdf(paste0(out.path,"/RegionRemoval.pdf"), width = 10, height = 9)
    graphics::par(mfrow=c(2,1), mar=c(4,4,2,2))
    xat=round(as.numeric(colnames(Spectrum_dataB)),5)
    graphics::plot(Re(Spectrum_dataB[nspectr,]), col="blue",xaxt="n", ylab = "Intensity", type="l", xlab="ppm", main="Before Region Removal")
    at=seq(1,length(xat), 30)
    graphics::axis(side=1, at=at, labels=round(xat[at],2))
    
    graphics::plot(Re(Spectrum_data[nspectr,]), col="blue",xaxt="n", ylab = "Intensity", type="l", xlab="ppm", main="After Region Removal")
    at=seq(1,length(xat), 30)
    graphics::axis(side=1, at=at, labels=round(xat[at],2))
#     text(x=441, y=160, labels=c("Lactate region"), srt = 90, col="red", cex=0.8)
#     text(x=269, y=100, labels=c("Water region"), col="red", cex=0.8)
#     arrows(x=255, y=70, x1 = 282, y1 = 70, length = 0.1,code = 3, col="red" )
#     arrows(x=439, y=70, x1 = 445, y1 = 70, length = 0.04,code = 3, col="red" )
    grDevices::dev.off()
  }
  
}


##########################
# ZoneAggregation
##########################

if (Za ==  TRUE ){
  Spectrum_dataB = Spectrum_data
  Spectrum_data = ZoneAggregation(Spectrum_data,fromto.za = fromto.za) 


  if (saveall == TRUE) {
    PretreatedSpectrabyStep[[step]] = Spectrum_data
    names(PretreatedSpectrabyStep)[step] <- "ZoneAggregation.Data"
    step = step + 1
  }
  
  if (ImpG==TRUE) {
    
    grDevices::pdf(paste0(out.path,"/ZoneAggregation.pdf"), width = 10, height = 9)
    graphics::par(mfrow=c(2,1), mar=c(4,4,2,2) )
    xat=round(as.numeric(colnames(Spectrum_data12)),5)
    graphics::plot(Re(Spectrum_dataB[nspectr,]), col="blue",  xaxt="n", ylab = "Intensity", type="l", xlab="ppm", main="Before Region Removal")
    at=seq(1,length(xat), 30)
    graphics::axis(side=1, at=at, labels=round(xat[at],2))
#     text(x=377, y=180, labels=c("Citrate region"),  col="red", cex=0.8)
#     arrows(x=373, y=150, x1 = 383, y1 = 150, length = 0.04,code = 3, col="red" )
    xat=round(as.numeric(colnames(Spectrum_data)),5)
    graphics::plot(Re(Spectrum_data[nspectr,]), col="blue", xaxt="n", ylab = "Intensity", type="l", xlab="ppm", main="Before Region Removal")
    at=seq(1,length(xat), 30)
    graphics::axis(side=1, at=at, labels=round(xat[at],2))
    grDevices::dev.off()
  }
  
}


##########################
# Normalization
##########################
if (N ==  TRUE ){

  Spectrum_dataB = Spectrum_data
  Spectrum_data = Normalization(Spectrum_data, type.norm=type.norm, from.norm=from.norm, to.norm=to.norm, ref.norm = ref.norm)


  if (saveall == TRUE) {
    PretreatedSpectrabyStep[[step]] = Spectrum_data
    names(PretreatedSpectrabyStep)[step] <- "Normalization.Data"
    step = step + 1
  }
  
  
  if (ImpG==TRUE) {
    grDevices::pdf(paste0(out.path,"/Normalization1.pdf"), width = 10, height = 5)
    graphics::par(mfrow=c(2,1))
    graphics::plot(Re(Spectrum_dataB[nspectr,]), col="blue", xaxt="n", type="l", xlab="ppm",  main="Before Normalization \n Real part")
    xat=round(as.numeric(colnames(Spectrum_data)),5)
    at=seq(1,length(xat), length(xat)/20)
    graphics::axis(side=1, at=at, labels=round(xat[at],2))
    graphics::lines(Re(Spectrum_data[nspectr,]), col= "red")
    graphics::legend("topleft", lty = 1, legend = c("Before Normalization", "After Normalization"), col = c("black", "red"))
    

    grDevices::pdf(paste0(out.path,"/Normalization2.pdf"), width = 10, height = 9)
    
    if (dim(Spectrum_dataB)[1]>=3) {
      f=sample(c(1:dim(Spectrum_dataB)[1]), 3)
    }else {f = 1}
    
    aa = c(0.5)*length(Spectrum_dataB[nspectr,])
    bb = c(0.75)*length(Spectrum_dataB[nspectr,])
    graphics::par(mfrow=c(2,1), mar=c(4,4,2,2))
    
    xat=round(as.numeric(colnames(Spectrum_dataB[,aa:bb])),5)
    graphics::plot(Re(Spectrum_dataB[nspectr,aa:b]), col=grDevices::rainbow(5,start=0.4)[1], xaxt="n", type="l",  ylab = "Intensity",xlab="ppm",  main="Before Normalization (zoom)")
    at=seq(1,length(xat), length(xat)/20)
    graphics::axis(side=1, at=at, labels=round(xat[at],2))
    for (i in 1:4){
      graphics::lines(Re(Spectrum_dataB[f[i],aa:bb]), col=grDevices::rainbow(5,start=0.4)[i+1]) 
    }
    xat=round(as.numeric(colnames(Spectrum_data[,aa:b])),5)
    graphics::plot(Re(Spectrum_data[nspectr,aa:bb]), col=grDevices::rainbow(5,start=0.4)[1], xaxt="n", type="l",  ylab = "Intensity",xlab="ppm",  main="After Normalization (zoom)")
    at=seq(1,length(xat), length(xat)/20)
    graphics::axis(side=1, at=at, labels=round(xat[at],2))
    for (i in 1:4){
      graphics::lines(Re(Spectrum_data[f[i],aa:bb]), col=grDevices::rainbow(5,start=0.4)[i+1]) 
    }
    grDevices::dev.off()
  }
  
  
}



#################### 
# save THE RESULTS
####################

Arguments <- c(as.list(environment()))

# remove unwanted variables in the arguments list
Arguments[c("Fid_info","Fid_dataB","PretreatedSpectrabyStep","RawSpect_dataB","Fid_data","FourierTransformData",
"SolventRe","baseline","ApodFactor","Ss.res" ,"Apod.res","aa","bb","cc","at","xat")] = NULL



if (save == TRUE) {
  if (RetArgs == TRUE) {
    save(Spectrum_data, Arguments, Fid_info, file=paste0(out.path, "/",dataname , "_FinalResPC.RData"))
    } else {save(Spectrum_data, Fid_info, file=paste0(out.path, "/",dataname , "_FinalSpectra.RData"))}
}

 
if (saveall == TRUE) {
  if (RetArgs == TRUE) {
    save(PretreatedSpectrabyStep, Arguments, Fid_info, file=paste0(out.path, "/",dataname , "PretreatedSpectrabyStep.RData"))
    } else {save(PretreatedSpectrabyStep, Fid_info, file=paste0(out.path, "/",dataname , "PretreatedSpectrabyStep.RData"))}
}


if (export == "csv") {
  write.csv(Spectrum_data, file = paste0(out.path, "/",dataname , "_Spectrum_data.csv"))
  write.csv(Fid_info, file = paste0(out.path, "/",dataname , "_FidInfo.csv"))
  }else if (export == "rdata") {
    if (RetArgs == TRUE) {save(Spectrum_data, Arguments, Fid_info, file = paste0(out.path, "/",dataname , "_FinalResPC.RData"))}
    save(Spectrum_data, Fid_info, file = paste0(out.path, "/",dataname , "_FinalResPC.RData"))
  }


if (RetArgs == TRUE) {
  return(list(Spectrum_data = Spectrum_data, Arguments = Arguments, Fid_info = Fid_info))
} else {return(list(Spectrum_data = Spectrum_data, Fid_info = Fid_info))}


}


