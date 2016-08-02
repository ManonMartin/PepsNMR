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
                        nspectr = 1, save = FALSE, saveall = FALSE, ImpG= FALSE, RetArgs = TRUE,
                        Fopc = TRUE, Ss = TRUE, A = TRUE, Zopc = TRUE, Bc = TRUE, 
                        Zsnv = TRUE, W = TRUE, B = TRUE, Zs = TRUE, Za=FALSE, N = TRUE,
                        l=1, subdirs = FALSE, #ReadFids
                        group_delay=NULL, # FOPC
                        lambda.ss=1e6, ptw.ss=TRUE,  # SolventSuppression
                        DT=NULL,type.apod = "exp",phase=0, rectRatio=1/2, gaussLB=1, expLB=1, # Apodization
                        SW_h=NULL, # FourierTransform
                        ptw.bc=TRUE, maxIter = 42,lambda.bc=1e7, p=0.05, eps=1e-8, # BaselineCorrection
                        shiftHandling="cut", thres = 300, #PPMConversion
                        
                        normalization.type="median", from.normW=3.05, to.normW=4.05,reference.choosing="fixed", 
                        reference=1,optim.crit="RMS", ptw.wp=F, K=3, L=40,
                        lambda.smooth=0, deg=3, lambda.bspline=0.01, kappa=0.0001,
                        max_it_Bspline=10, returnReference=F,  #Warping
                        from.ws = 0.2, to.ws = 10, reverse.axis = TRUE, # WindowSelection
                        m = 500, # Bucketing
                        typeofspectra = "serum",type.rr =  "zero", fromto.rr=list(Water =c(4.5, 5.1)), # RegionRemoval
                        fromto.za = list(Citrate =c(2.5, 2.7)), # ZoneAggregation
                        type.norm="mean", from.norm=3.05, to.norm=4.05, ref.norm # Normalization
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
  longueur = length(which(c(Fopc , Ss , A , Zopc, Bc ,Zsnv , W , B, Zs , Za, N)==TRUE)) # optional steps
  longueur = longueur + 4 # nombre d'etapes totales
    
  PretreatedSpectrabyStep = vector("list", longueur) # creation d une liste pour sauver les donnees de chaque etape
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

  Fid_data1 = Fid_data
  
  
  if (saveall == TRUE) {
    PretreatedSpectrabyStep[[step]] = Fid_data1 
    names(PretreatedSpectrabyStep)[step] <- "FirstOrderPhaseCorrection.Data"
    step = step + 1
  }
 

  if (ImpG==TRUE) {
    grDevices::pdf(paste0(out.path,"/FirstOrderPhCorrectREAL.pdf"),width=13, height=8)
    graphics::par(mfrow=c(2,2))
    graphics::plot(Re(Fid_dataB[nspectr,]), col="blue", type="l", ylab = "Intensity", xlab="Time", main="Original FID \n Real part")
    graphics::plot(Re(Fid_dataB[nspectr,1:200]), col="blue", ylab = "Intensity",xlab="Time", type="l", main="Original FID (zoom)\n Real part ")
    graphics::plot(Re(Fid_data1[nspectr,]), col="blue",  ylab = "Intensity",xlab="Time", type="l", main="1st Order Phase Correction \n Real part")
    graphics::plot(Re(Fid_data1[nspectr,1:200]), col="blue", ylab = "Intensity",xlab="Time", type="l", main="1st Order Phase Correction (zoom)\n Real part")
    grDevices::dev.off()
    
    grDevices::pdf(paste0(out.path,"/FirstOrderPhCorrectIMAG.pdf"),width=13, height=8)
    graphics::par(mfrow=c(2,2))
    graphics::plot(Im(Fid_dataB[nspectr,]), col="blue", ylab = "Intensity",xlab="Time", type="l", main="Original FID \n Imaginary part")
    graphics::plot(Im(Fid_dataB[nspectr,1:200]), col="blue", ylab = "Intensity",xlab="Time", type="l", main="Original FID (zoom)\n Imaginary part ")
    graphics::plot(Im(Fid_data1[nspectr,]), col="blue", ylab = "Intensity",xlab="Time", type="l", main="1st Order Phase Correction \n Imaginary part")
    graphics::plot(Im(Fid_data1[nspectr,1:200]), col="blue", ylab = "Intensity",xlab="Time", type="l", main="1st Order Phase Correction (zoom)\n Imaginary part ")
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
  
  Fid_data2 = Fid_data

  
  if (saveall == TRUE) {
    PretreatedSpectrabyStep[[step]] = Fid_data2 
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
    graphics::plot(Re(Fid_data2[nspectr,1:100]), col="blue", type="l",  ylab = "Intensity", xlab="Time", main="After solvent suppression (zoom)\n Real part ")
    graphics::plot(Im(Fid_data2[nspectr,1:100]), col="blue", type="l",  ylab = "Intensity", xlab="Time", main="After solvent suppression (zoom)\n Imaginary part ")
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
  
  Fid_data3 = Fid_data


  if (saveall == TRUE) {
    PretreatedSpectrabyStep[[step]] = Fid_data3 
    names(PretreatedSpectrabyStep)[step] <- "Apodization.Data"
    step = step + 1
  }
  
  if (ImpG==TRUE) {
    grDevices::pdf(paste0(out.path,"/Apodization1.pdf"),width=13, height=8)
    graphics::par(mfrow=c(2,2))
    graphics::plot(Re(Fid_dataB[nspectr,]), col="blue", type="l", xlab="Time", main="Before apodization \n Real part")
    graphics::plot(Im(Fid_dataB[nspectr,]), col="blue", type="l", xlab="Time", main="Before apodization \n Imaginary part")
    graphics::plot(Re(Fid_data3[nspectr,]), col="blue", type="l", xlab="Time", main="After apodization \n Real part")
    graphics::plot(Im(Fid_data3[nspectr,]), col="blue", type="l", xlab="Time", main="After apodization \n Imaginary part")
    grDevices::dev.off()
    
    grDevices::pdf(paste0(out.path,"/Apodization2.pdf"),width=8, height=4)
    graphics::plot(Apod.res[["Factor"]], col="red",type="l", main = "Apodization factor")
    grDevices::dev.off()
  }

}

##########################
# FourierTransform
##########################

Fid_dataB = Fid_data
RawSpect_data = FourierTransform(Fid_data, Fid_info = Fid_info, SW_h = SW_h)

RawSpect_data4 = RawSpect_data
FourierTransformData = Re(RawSpect_data)

if (save == TRUE) {
  save(FourierTransformData = FourierTransformData,  file=paste0(out.path, "/",dataname ,"FourierTransformData.RData"))
}

if (saveall == TRUE) {
  PretreatedSpectrabyStep[[step]] = RawSpect_data4 
  names(PretreatedSpectrabyStep)[step] <- "FourierTransform.Data"
  step = step + 1
}


if (ImpG==TRUE) {
  grDevices::pdf(paste0(out.path,"/FourierTransform.pdf"),width=13, height=14)
  graphics::par(mfrow=c(3,1))
  a = 0.2*length(RawSpect_data4[nspectr,])
  b = 0.6*length(RawSpect_data4[nspectr,])
  c = 0.06*length(RawSpect_data4[nspectr,])
  
  graphics::plot(Re(RawSpect_data4[nspectr,]), col="blue", type="l", xlab="Frequency", main="After Fourier Transform \n Real part")
  graphics::plot(Re(RawSpect_data4[nspectr,a:b]), col="blue", xaxt = "n", type="l", xlab="Frequency", main="After Fourier Transform \n Real part (zoom)")
  at=seq(1,(b-a), c)
  graphics::axis(side=1, at=at, labels  = c(a:b)[at])
  graphics::plot(Im(RawSpect_data4[nspectr,a:b]), col="blue", xaxt = "n", type="l", xlab="Frequency", main="After Fourier Transform \n Imaginary part (zoom)")
  at=seq(1,(b-a), c)
  graphics::axis(side=1, at=at, labels = c(a:b)[at])
  grDevices::dev.off()
  ###
}

##########################
# ZeroOrderPhaseCorrection
##########################


if (Zopc ==  TRUE ){

  RawSpect_dataB = RawSpect_data
  RawSpect_data = ZeroOrderPhaseCorrection(RawSpect_data, plot_rms = NULL)

  RawSpect_data5 = RawSpect_data


  if (saveall == TRUE) {
    PretreatedSpectrabyStep[[step]] = RawSpect_data5 
    names(PretreatedSpectrabyStep)[step] <- "ZeroOrderPhaseCorrection.Data"
    step = step + 1
  }
  
  
  if (ImpG==TRUE) {
    grDevices::pdf(paste0(out.path,"/ZeroOrderPhaseCorrection.pdf"),width=13, height=14)
    graphics::par(mfrow=c(3,1))
    
    a = 0.2*length(RawSpect_dataB[nspectr,])
    b = 0.6*length(RawSpect_dataB[nspectr,])
    c = 0.06*length(RawSpect_dataB[nspectr,])
    
    graphics::plot(Re(RawSpect_data4[nspectr,a:b]), col="blue",  xaxt = "n", type="l", xlab="Frequency", main="Before 0 Order Phase Correction \n Real part (zoom)")
    at=seq(1,(b-a), c)
    graphics::axis(side=1, at=at, labels = c(a:b)[at])
    graphics::plot(Re(RawSpect_data5[nspectr,a:b]), col="blue",  xaxt = "n", type="l", xlab="Frequency", main="After 0 Order Phase Correction \n Real part (zoom)")
    at=seq(1,(b-a), c)
    graphics::axis(side=1, at=at, labels = c(a:b)[at])
    graphics::plot(Im(RawSpect_data5[nspectr,a:b]), col="blue",  xaxt = "n", type="l", xlab="Frequency", main="After 0 Order Phase Correction \n Imaginary part (zoom)")
    at=seq(1,(b-a), c)
    graphics::axis(side=1, at=at, labels = c(a:b)[at])
    grDevices::dev.off()
    }

}

##########################
# BaselineCorrection
##########################
if (Bc ==  TRUE ){

  RawSpect_dataB = RawSpect_data
  BC.res =  BaselineCorrection(RawSpect_data, returnBaseline=TRUE, lambda.bc=lambda.bc, p=p, eps=eps)
  baseline = BC.res[["Baseline"]]
  RawSpect_data = BC.res[["RawSpect_data"]]
  
  RawSpect_data6 = RawSpect_data



  if (saveall == TRUE) {
    PretreatedSpectrabyStep[[step]] = RawSpect_data6 
    names(PretreatedSpectrabyStep)[step] <- "BaselineCorrection.Data"
    step = step + 1
  }

  if (ImpG==TRUE) {
    grDevices::pdf(paste0(out.path,"/BaselineCorrection.pdf"),width=13, height=8)
    graphics::par(mfrow=c(2,2))
    
    a = 0.42*length(RawSpect_dataB[nspectr,])
    b = 0.46*length(RawSpect_dataB[nspectr,])
    
    graphics::plot(Re(RawSpect_dataB[nspectr,]), col="blue", type="l", ylab = "Intensity", xlab="Frequency", main="Before Baseline Correction")
    graphics::abline(h=0, lty = 2)
    graphics::lines(baseline[nspectr,], col="red")
    
    graphics::plot(Re(RawSpect_dataB[nspectr,a:b]), col="blue", xaxt="n", type="l", ylab = "Intensity", xlab="Frequency", main="Before Baseline Correction  \n zoom")
    graphics::axis(side = 1, at = seq(0,(b-a), 100), labels = seq(a,b, 100))
    graphics::abline(h=0, lty = 2)
    graphics::lines(baseline[nspectr,a:b], col="red")
    
    graphics::plot(Re(RawSpect_data6[nspectr,]), col="blue", type="l", ylab = "Intensity", xlab="Frequency", main="After BaselineCorrection")
    graphics::abline(h=0, lty = 2)
    
    graphics::plot(Re(RawSpect_data6[nspectr,a:b]),  col="blue", xaxt="n", type="l", ylab = "Intensity", xlab="Frequency", main="After BaselineCorrection \n zoom")
    graphics::axis(side = 1, at = seq(0,(b-a), 100), labels = seq(a,b, 100))
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
  
  RawSpect_data7 = RawSpect_data


  if (saveall == TRUE) {
    PretreatedSpectrabyStep[[step]] = RawSpect_data7 
    names(PretreatedSpectrabyStep)[step] <- "NegativeValuesZeroing.Data"
    step = step + 1
  }
  
  if (ImpG==TRUE) {
    grDevices::pdf(paste0(out.path,"/NegativeValuesZeroing.pdf"), width = 13, height = 8)
    graphics::par(mfrow=c(2,1))
    
    a = 0.2*length(RawSpect_dataB[nspectr,])
    b = 0.6*length(RawSpect_dataB[nspectr,])
    c = 0.06*length(RawSpect_dataB[nspectr,])
    
    graphics::plot(Re(RawSpect_dataB[nspectr,a:b]), col="blue", xaxt = "n",type="l", xlab="Frequency", main="Before Negative Values Zeroing \n Real part (zoom)")
    at=seq(1,(b-a), c)
    graphics::axis(side=1, at=at, labels = c(a:b)[at])
    graphics::plot(Re(RawSpect_data7[nspectr,a:b]), col="blue", xaxt = "n",type="l", xlab="Frequency", main="After Negative Values Zeroing \n Real part (zoom)")
    at=seq(1,(b-a), c)
    graphics::axis(side=1, at=at, labels = c(a:b)[at])
    grDevices::dev.off()
  }
  
}


##########################
# PPMConversion 
##########################
RawSpect_dataB = RawSpect_data
Spectrum_data = PPMConversion(RawSpect_data, RawSpect_info = Fid_info, shiftHandling=shiftHandling, thres = thres)

Spectrum_data8 = Spectrum_data


if (saveall == TRUE) {
  PretreatedSpectrabyStep[[step]] = Spectrum_data8 
  names(PretreatedSpectrabyStep)[step] <- "PPMConversion.Data"
  step = step + 1
}

if (ImpG==TRUE) {
  grDevices::pdf(paste0(out.path,"/PPMConversion.pdf"),width=13, height=8)
  graphics::par(mfrow=c(2,1))
  graphics::plot(Re(RawSpect_dataB[nspectr,]), col="blue", type="l", xlab="Frequency", main="Before PPM Conversion \n Real part")
  xat=round(as.numeric(colnames(Spectrum_data8)),5)
  graphics::plot(xat,Re(Spectrum_data8[nspectr,]), col="blue", type="l", xlab="ppm", main="After PPM Conversion \n Real part")
  grDevices::dev.off()
}


############
# Warping
############

if (W ==  TRUE ){
  Spectrum_dataB = Spectrum_data
  Warp.res  = Warping(RawSpect_data=Spectrum_data, returnWarpingfunc=TRUE, normalization.type=normalization.type, 
                      from.normW=from.normW, to.normW=to.normW,reference.choosing=reference.choosing, 
                      reference=reference,optim.crit=optim.crit, ptw.wp=ptw.wp, K=K, L=L,
                      lambda.smooth=lambda.smooth, deg=deg, lambda.bspline=lambda.bspline, kappa=kappa,
                      max_it_Bspline=max_it_Bspline, returnReference=returnReference)
  Spectrum_data =  Warp.res[["RawSpect_data"]]
  warpingfunc = Warp.res[["Warpingfunc"]]

  Spectrum_data9 = Spectrum_data
  
  
  if (saveall == TRUE) {
    PretreatedSpectrabyStep[[step]] = Spectrum_data9 
    names(PretreatedSpectrabyStep)[step] <- "Warping.Data"
    step = step + 1
  }
  
  if (ImpG==TRUE) {
    grDevices::pdf(paste0(out.path,"/Warping1.pdf"), width = 8, height = 4)
    graphics::plot(warpingfunc, col="blue", type="l", main = "warping function")
    grDevices::dev.off()
    
    grDevices::pdf(paste0(out.path,"/Warping2.pdf"), width = 13, height = 10)
    graphics::par(mfrow=c(2,1),mar=c(4.1, 4.1, 4.1, 9.5), xpd=TRUE)
    
    if (dim(Spectrum_dataB)[1]>=3) {
      f=sample(c(1:dim(Spectrum_dataB)[1]), 3)
    }else {f = 1}

    
    a = 0.4*length(RawSpect_dataB[f[1],])
    b = 0.47*length(RawSpect_dataB[f[1],])
    c = 0.015*length(RawSpect_dataB[f[1],])
    
    graphics::plot(Re(Spectrum_dataB[nspectr,a:b]),  col="red", xaxt = "n",ylab = "Intensity",ylim=c(0, max(Re(Spectrum_dataB[c(nspectr, f),a:b]))), type="l", xlab="Frequency", main="Before Warping (zoom)")
    graphics::axis(side = 1, at = seq(0,(b-a), c), labels = seq(a,b, c))   
    graphics::legend("topright", inset=c(-0.22,-0.12), legend=c("Ref. spectrum", "Warped  spectra"), lty = 1, cex=0.8, col=c("red", "blue"))    
    for (j in f) {
      graphics::lines(Re(Spectrum_dataB[j,a:b]), col="blue", type="l")}
    graphics::grid(20, NA, lty = 1, lwd = 0.7, col="gray44")
    
    graphics::plot(Re(Spectrum_data9[nspectr,a:b]), col="red", xaxt = "n",ylab = "Intensity",ylim=c(0, max(Re(Spectrum_data9[c(nspectr, f),a:b]))), type="l", xlab="Frequency", main="After Warping (zoom)")
    graphics::axis(side = 1, at = seq(0,(b-a), c), labels = seq(a,b, c))    
    graphics::legend("topright", inset=c(-0.22,-0.12), legend=c("Ref. spectrum", "Warped  spectra"), lty = 1, cex=0.8, col=c("red", "blue"))    
    for (j in f) {
      graphics::lines(Re(Spectrum_data9[j,a:b]), type="l", col="blue")}
    graphics::grid(20, NA, lty = 1, lwd = 0.7, col="gray44")
    grDevices::dev.off()
  }
  
}


##########################
# WindowSelection
##########################
Spectrum_dataB = Spectrum_data
Spectrum_data = WindowSelection(Spectrum_data, from.ws = from.ws, to.ws = to.ws, reverse.axis = reverse.axis)

Spectrum_data10 = Spectrum_data

if (saveall == TRUE) {
  PretreatedSpectrabyStep[[step]] = Spectrum_data10 
  names(PretreatedSpectrabyStep)[step] <- "WindowSelection.Data"
  step = step + 1
}


if (ImpG==TRUE) {
  grDevices::pdf(paste0(out.path,"/WindowSelection.pdf"), width = 10, height = 9)
  graphics::par(mfrow=c(2,1), mar=c(4,4,2,2))
  c = 0.03*length(RawSpect_dataB[nspectr,])
  xat=round(as.numeric(colnames(Spectrum_dataB)),5)
  graphics::plot(xat,Re(Spectrum_dataB[nspectr,]),col="blue",  type="l",  ylab = "Intensity", xlab="ppm", main="Before Window Selection")
  xat=round(as.numeric(colnames(Spectrum_data10)),5)
  graphics::plot( Re(Spectrum_data10[nspectr,]), col="blue",  ylab = "Intensity", xaxt="n", type="l", xlab="ppm", main="After Window Selection")
  at=seq(1,length(xat), c)
  graphics::axis(side=1, at=at, labels=round(xat[at],2))
  grDevices::dev.off()
}
 
################ 
# Bucketing
################

if (B ==  TRUE ){
  Spectrum_dataB = Spectrum_data
  Spectrum_data = Bucketing(Spectrum_data, m = m)
  
  Spectrum_data11 = Spectrum_data

  if (saveall == TRUE) {
    PretreatedSpectrabyStep[[step]] = Spectrum_data11 
    names(PretreatedSpectrabyStep)[step] <- "Bucketing.Data"
    step = step + 1
  }
  
  if (ImpG==TRUE) {
    grDevices::pdf(paste0(out.path,"/Bucketing.pdf"), width = 10, height = 9)
    graphics::par(mfrow=c(2,1), mar=c(4,4,2,2))
    c = 0.03*length(RawSpect_dataB[nspectr,])
    xat=round(as.numeric(colnames(Spectrum_dataB)),5)
    graphics::plot(Re(Spectrum_dataB[nspectr,]), col="blue",ylab = "Intensity", xaxt="n", type="l", xlab="ppm", main="Before Bucketing")
    at=seq(1,length(xat), c)
    graphics::axis(side=1, at=at, labels=round(xat[at],2))
    xat=round(as.numeric(colnames(Spectrum_data11)),5)
    graphics::plot(Re(Spectrum_data11[nspectr,]), col="blue",ylab = "Intensity", xaxt="n", type="l", xlab="ppm", main="After Bucketing")
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

  Spectrum_data12 = Spectrum_data
  

  if (saveall == TRUE) {
    PretreatedSpectrabyStep[[step]] = Spectrum_data12 
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
    
    graphics::plot(Re(Spectrum_data12[nspectr,]), col="blue",xaxt="n", ylab = "Intensity", type="l", xlab="ppm", main="After Region Removal")
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

  
  Spectrum_data13 = Spectrum_data


  if (saveall == TRUE) {
    PretreatedSpectrabyStep[[step]] = Spectrum_data13 
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
    xat=round(as.numeric(colnames(Spectrum_data13)),5)
    graphics::plot(Re(Spectrum_data13[nspectr,]), col="blue", xaxt="n", ylab = "Intensity", type="l", xlab="ppm", main="Before Region Removal")
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
  
  Spectrum_data14 = Spectrum_data



  if (saveall == TRUE) {
    PretreatedSpectrabyStep[[step]] = Spectrum_data14 
    names(PretreatedSpectrabyStep)[step] <- "Normalization.Data"
    step = step + 1
  }
  
  
  if (ImpG==TRUE) {
    grDevices::pdf(paste0(out.path,"/Normalization1.pdf"), width = 10, height = 5)
    graphics::par(mfrow=c(2,1))
    graphics::plot(Re(Spectrum_dataB[nspectr,]), col="blue", xaxt="n", type="l", xlab="ppm",  main="Before Normalization \n Real part")
    xat=round(as.numeric(colnames(Spectrum_data14)),5)
    at=seq(1,length(xat), length(xat)/20)
    graphics::axis(side=1, at=at, labels=round(xat[at],2))
    graphics::lines(Re(Spectrum_data14[nspectr,]), col= "red")
    graphics::legend("topleft", lty = 1, legend = c("Before Normalization", "After Normalization"), col = c("black", "red"))
    

    grDevices::pdf(paste0(out.path,"/Normalization2.pdf"), width = 10, height = 9)
    
    if (dim(Spectrum_dataB)[1]>=3) {
      f=sample(c(1:dim(Spectrum_dataB)[1]), 3)
    }else {f = 1}
    
    a = c(0.5)*length(Spectrum_dataB[nspectr,])
    b = c(0.75)*length(Spectrum_dataB[nspectr,])
    graphics::par(mfrow=c(2,1), mar=c(4,4,2,2))
    
    xat=round(as.numeric(colnames(Spectrum_dataB[,a:b])),5)
    graphics::plot(Re(Spectrum_dataB[nspectr,a:b]), col=grDevices::rainbow(5,start=0.4)[1], xaxt="n", type="l",  ylab = "Intensity",xlab="ppm",  main="Before Normalization (zoom)")
    at=seq(1,length(xat), length(xat)/20)
    graphics::axis(side=1, at=at, labels=round(xat[at],2))
    for (i in 1:4){
      graphics::lines(Re(Spectrum_dataB[f[i],a:b]), col=grDevices::rainbow(5,start=0.4)[i+1]) 
    }
    xat=round(as.numeric(colnames(Spectrum_data14[,a:b])),5)
    graphics::plot(Re(Spectrum_data14[nspectr,a:b]), col=grDevices::rainbow(5,start=0.4)[1], xaxt="n", type="l",  ylab = "Intensity",xlab="ppm",  main="After Normalization (zoom)")
    at=seq(1,length(xat), length(xat)/20)
    graphics::axis(side=1, at=at, labels=round(xat[at],2))
    for (i in 1:4){
      graphics::lines(Re(Spectrum_data14[f[i],a:b]), col=grDevices::rainbow(5,start=0.4)[i+1]) 
    }
    grDevices::dev.off()
  }
  
  
}



#################### 
# save THE RESULTS
####################
argnames <- c("dataname",           "data.path",          "out.path", 
              "nspectr",            "save",           
              "saveall",            "ImpG",               "Fopc",              
              "Ss" ,                "A",                  "Zopc",               "Bc",                
              "Zsnv" ,              "W",                  "B",                  "Zs" ,               
              "Za",                 "N",                   "l",        "subdirs" ,          
              "group_delay",        "lambda.ss",          "ptw.ss",                   
              "DT",                 "type.apod",          "phase",              "rectRatio",         
              "gaussLB",            "expLB",              "SW_h",               "",          
              "ptw.bc",             "maxIter",            "lambda.bc",         
              "p",                  "eps",     "",           "shiftHandling",      "thres",                 
              "normalization.type", "from.normW",         "to.normW",          
              "reference.choosing", "reference",          "optim.crit",         "ptw.wp",            
              "K",                  "L",                  "lambda.smooth",      "deg"  ,             
              "lambda.bspline",     "kappa",              "max_it_Bspline",     "returnReference" ,  
              "from.ws",            "to.ws",              "reverse.axis",       "m"  ,               
              "typeofspectra",      "type.rr",            "fromto.rr",          "fromto.za",         
              "type.norm",          "from.norm",          "to.norm",            "ref.norm")

supargnames = c(rep("general", 18), 
                rep("ReadFids", 2), 
                rep("FirstOrderPhaseCorrection", 1), 
                rep("SolventSuppression",2), 
                rep("Apodization",6), 
                rep("FourierTransform", 1),
                
                rep("Zero order phase correction", 1),
                rep("BaselineCorrection", 5),
                
                rep("Negative values zeroing", 1),
                rep("PPMConversion", 2),
                rep("Warping", 15),
                rep("WindowSelection", 3),
                rep("Bucketing", 1),
                rep("RegionRemoval", 3),
                rep("ZoneAggregation", 1),
                rep("Normalization", 4))



Spectra=Re(Spectrum_data)

Arguments <- c(as.list(environment()))
index <- names(Arguments) %in% argnames
Arguments <- Arguments[index]

Arguments = data.frame(cbind(Category = supargnames, Argument_name = argnames, Argument_value = c(Arguments[1:31], "",Arguments[32:36], "" ,Arguments[37:66])), row.names = NULL)



if (save == TRUE) {
  if (RetArgs == TRUE) {
    save(Spectra, Arguments, Fid_info, file=paste0(out.path, "/",dataname , "FinalRes.RData"))
    } else {save(Spectra, file=paste0(out.path, "/",dataname , "FinalSpectra.RData"))}
}
  
if (saveall == TRUE) {
save(PretreatedSpectrabyStep, file=paste0(out.path, "/",dataname , "PretreatedSpectrabyStep.RData"))
}


if (RetArgs == TRUE) {
  return(list(Spectra = Spectra, Arguments = Arguments, Fid_info = Fid_info))
} else {return(Spectra)}


}
