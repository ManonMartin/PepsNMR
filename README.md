[![Build Status](https://travis-ci.org/ManonMartin/PepsNMR.svg?branch=master)](https://travis-ci.org/ManonMartin/PepsNMR)

# PepsNMR
An R package for 1H NMR data pre-processing.

## R code to install and load/attach the package:

#### First option:

```R
require(devtools)
install_github("ManonMartin/PepsNMR", dependencies = TRUE)
require(PepsNMR)
```

#### Second option:

Download first the compressed package files from GithHub, unzip it and save the local copy, then:

```R
install.packages("path-to-PepsNMR_folder", repos = NULL, type="source")
require(PepsNMR)
```

*Rem*: if you use the second option, ensure that all the dependencies (i.e. the packages listed in the DESCRIPTION (Imports:) and NAMESPACE files (import and importFrom)) are correclty installed on your computer.


## Minimal example

The following R code is the application of the pre-processing steps to the Human Serum dataset included in PepsNMR with graphical outputs. It can be copy-pasted to pre-process your own data. 

```R
###########################################
#         A minimal example               #
###########################################

# ==== load/attach the package =================
require(PepsNMR)

# ==== set graphical parameters =================
# save default graphical parameters
default.par <- par(no.readonly=T)
# select the index of the spectrum that will be drawn
spectrIndex <- 1
# colors
col1 <-  "gray18"
col2 <- "firebrick1"

# ==== define the path to the data files =================
data_path <-  system.file("extdata", package = "PepsNMR")

# ==== read the FIDs and their metadata =================
fidList <- ReadFids(file.path(data_path, "HumanSerum"))
Fid_data0 <- fidList[["Fid_data"]]
Fid_info <- fidList[["Fid_info"]]

time <- as.numeric(colnames(Fid_data0))*1000
plot(time, Re(Fid_data0[spectrIndex,]),type="l", col = col2, xlab=expression(paste("Time (", 10^3*mu,"s)")), 
     ylab = "Intensity", main = "Raw FID (real part)", ylim = c(-1e6,7e5))


# ==== GroupDelayCorrection =================
Fid_data.GDC <- GroupDelayCorrection(Fid_data0, Fid_info)

par(mfrow=c(2,1), mar = c(4,4,2,2))
plot(time[0:300], Re(Fid_data0[spectrIndex,0:300]),  
     type="l", ylab = "Intensity", xlab="", 
     main="FID with the Group Delay (real part - zoom)", col = col1)
plot(time[0:300], Re(Fid_data.GDC[spectrIndex,0:300]), 
     type="l", ylab = "Intensity", xlab=expression(paste("Time (", mu,"s)")), 
     main="FID without the Group Delay (real part - zoom)", col = col1)


# ====  SolventSuppression =================
SS.res <- SolventSuppression(Fid_data.GDC, returnSolvent=TRUE)
Fid_data.SS <- SS.res[["Fid_data"]]
SolventRe <- SS.res[["SolventRe"]]

par(mar=c(4,4,1.5,1), mfrow=c(2,1))
plot(time[0:4000], Re(Fid_data.GDC[spectrIndex,0:4000]),  col=col1, 
     type="l", ylab = "Intensity", xlab="", 
     main="FID and solvent residuals signal (real part - zoom)")
lines(time[0:4000],SolventRe[spectrIndex,0:4000], col=col2 , lwd = 1.3)
legend("topright", bty = "n", legend = c(expression(paste("Estimated solvent residuals signal ", (italic(W)))), expression(paste("FID signal ", (italic(S))))), 
       col=c(col2, col1),  lty = 1)

plot(time[0:4000], Re(Fid_data.SS[1,0:4000]), col=col1, 
     type="l", ylab = "Intensity", xlab=expression(paste("Time (", mu,"s)")), 
     main="FID without solvent residuals signal (real part - zoom)")
lines(time[0:4000], rep(0, 4000), col=col2)



# ==== Apodization =================
Fid_data.A <- Apodization(Fid_data.SS, Fid_info)

par(mar=c(4,4,1.5,1), mfrow=c(2,1))
plot(time, Re(Fid_data.SS[spectrIndex,]),  col=col1, 
     type="l", ylab = "Intensity", xlab="", main="FID before Apodisation")

plot(time, Re(Fid_data.A[spectrIndex,]), col=col1, 
     type="l", ylab = "Intensity", xlab=expression(paste("Time (", mu,"s)")), 
     main="FID after Apodisation")



# ==== FourierTransform =================
RawSpect_data.FT <- FourierTransform(Fid_data.A, Fid_info)

par(default.par) 
plot(Re(RawSpect_data.FT[spectrIndex,]), col=col1, xaxt="n",
     type="l", ylab = "Intensity", xlab = "ppm", 
     main="Spectrum after Fourier Transform")
at <- seq(1,dim(RawSpect_data.FT)[2], floor(dim(RawSpect_data.FT)[2]/10))
axis(side=1, at = at, 
     labels = round(as.numeric(colnames(RawSpect_data.FT)[at]),2))


# ==== ZeroOrderPhaseCorrection =================
Spectrum_data.ZOPC <- ZeroOrderPhaseCorrection(RawSpect_data.FT)

par(default.par) 
plot(Re(Spectrum_data.ZOPC[spectrIndex,]), col=col1, xaxt="n",
     type="l", ylab = "Intensity", xlab = "ppm", 
     main="Spectrum after Zero Order Phase Correction")
at <- seq(1,dim(Spectrum_data.ZOPC)[2], floor(dim(Spectrum_data.ZOPC)[2]/10))
axis(side=1, at = at, 
     labels = round(as.numeric(colnames(Spectrum_data.ZOPC)[at]),2))


# ==== InternalReferencing =================
target.value <- 0
Spectrum_data.IR <- InternalReferencing(Spectrum_data.ZOPC, Fid_info,
                                        ppm.value = target.value)

par(default.par) 
ppmvalues <- as.numeric(colnames(Spectrum_data.IR))
plot(Re(Spectrum_data.IR[spectrIndex,]), col=col1, xaxt="n",
     type="l", ylab = "Intensity", xlab = "ppm", 
     main="Spectrum after Internal Referencing")
at <- seq(1,dim(Spectrum_data.IR)[2], floor(dim(Spectrum_data.IR)[2]/10))
axis(side=1, at = at, 
     labels = round(ppmvalues[at],2))

index <- which(abs(ppmvalues-target.value) == min(abs(ppmvalues-target.value)))
abline(v = index, col= col2)
legend("topright", bty = "n", legend = "Peak location", 
       col=col2,  lty = 1)



# ==== BaselineCorrection =================
BC.res <- BaselineCorrection(Spectrum_data.IR, returnBaseline = TRUE,
                             lambda.bc = 1e8, p.bc = 0.01)

Spectrum_data.BC <- BC.res[["Spectrum_data"]] 
Baseline <- BC.res[["Baseline"]]

par(mar=c(4,4,1,1), mfrow=c(2,1))
plot(Re(Spectrum_data.IR[spectrIndex,]), col=col1, xaxt="n",
     type="l", ylab = "Intensity", xlab = "", 
     main="Spectrum before Baseline Correction")
at <- seq(1,dim(Spectrum_data.IR)[2], floor(dim(Spectrum_data.IR)[2]/10))
axis(side=1, at = at, labels = round(ppmvalues[at],2))
lines(Baseline[,1], col=col2)
legend("topright", bty = "n", legend = "Estimated baseline ", 
       col = col2,  lty = 1)

plot(Re(Spectrum_data.BC[spectrIndex,]), col=col1, xaxt="n",
     type="l", ylab = "Intensity", xlab = "ppm", 
     main="Spectrum after Baseline Correction")
axis(side=1, at = at, labels = round(ppmvalues[at],2))


# ==== NegativeValuesZeroing =================
Spectrum_data.NVZ <- NegativeValuesZeroing(Spectrum_data.BC)

par(default.par) 
plot(Re(Spectrum_data.NVZ[spectrIndex,]), col=col1, xaxt="n",
     type="l", ylab = "Intensity", xlab = "ppm", 
     main="Spectrum after Negative Values Zeroing")
axis(side=1, at = at, labels = round(ppmvalues[at],2))


# ==== Warping =================
W.res <- Warping(Spectrum_data.NVZ, returnWarpFunc = TRUE, 
                 reference.choice = "fixed")

Spectrum_data.W <- W.res[["Spectrum_data"]]
warp_func <- W.res[["Warp.func"]]

par(mfrow=c(2,1),mar=c(4,4,1.5,1))
f = c(21, 20, 24) # warped spectra index to draw
fen = c(17780:18240) # x-window
ylim = c(0, max(c(Re(Spectrum_data.NVZ[c(1, f),fen]), Re(Spectrum_data.W[c(spectrIndex, f),fen]))))

# Unwarped spectra
plot(Re(Spectrum_data.NVZ[1, fen]),   xaxt = "n", col=col2, ylab = "Intensity",ylim=ylim, type="l", xlab="ppm", main="Spectra before warping (real part - zoom)")
legend("topright", bty = "n", y.intersp = 0.8,legend=c("Unwarped spectra","Ref. spectrum "), lty = c(1,1), col=c(col1,col2))    
axis(side=1,  at = seq(1,length(fen), 114), labels = round(as.numeric(colnames(Spectrum_data.NVZ[,fen])[seq(1,length(fen), 114)]),2))
for (j in f) {
  graphics::lines(Re(Spectrum_data.NVZ[j,fen]), col=col1, type="l")
  }

# Warped spectra
plot(Re(Spectrum_data.W[1, fen]), col=col2, xaxt = "n",ylab = "Intensity",ylim=ylim, type="l", xlab="ppm", main="Warped spectra (real part - zoom)")
legend("topright",   bty = "n",  y.intersp = 0.8, legend=c("Warped spectra ","Ref. spectrum "), lty = c(1,1), col=c(col1,col2))    
axis(side=1,  at = seq(1,length(fen), 114), labels = round(as.numeric(colnames(Spectrum_data.NVZ[,fen])[seq(1,length(fen), 114)]),2))
for (j in f) {
  graphics::lines(Re(Spectrum_data.W[j,fen]), col=col1, type="l")
  }



# ==== WindowSelection =================
Spectrum_data.WS <- WindowSelection(Spectrum_data.W, from.ws = 10, to.ws = 0.2)

par(default.par) 
at <- seq(1,dim(Spectrum_data.WS)[2], floor(dim(Spectrum_data.WS)[2]/10))
ppmvalues <- as.numeric(colnames(Spectrum_data.WS))
plot(Re(Spectrum_data.WS[spectrIndex,]), col = col1, xaxt = "n",
     type = "l", ylab = "Intensity", xlab = "ppm", 
     main = "Spectrum after Window Selection")
axis(side = 1, at = at, labels = round(ppmvalues[at],2))



# ==== Bucketing =================
Spectrum_data.B <- Bucketing(Spectrum_data.WS, intmeth = "t")

par(mar=c(4,4,1,1), mfrow=c(2,1))
at <- seq(1,dim(Spectrum_data.WS)[2], floor(dim(Spectrum_data.WS)[2]/10))
ppmvalues <- as.numeric(colnames(Spectrum_data.WS))
plot(Re(Spectrum_data.WS[spectrIndex,]), col = col1, xaxt = "n",
     type = "l", ylab = "Intensity", xlab = "", 
     main = "Spectrum before Bucketing")
axis(side = 1, at = at, labels = round(ppmvalues[at],2))

at <- seq(1,dim(Spectrum_data.B)[2], floor(dim(Spectrum_data.B)[2]/10))
ppmvalues <- as.numeric(colnames(Spectrum_data.B))
plot(Re(Spectrum_data.B[spectrIndex,]), col = col1, xaxt = "n",
     type = "l", ylab = "Intensity", xlab = "ppm", 
     main = "Spectrum after Bucketing")
axis(side = 1, at = at, labels = round(ppmvalues[at],2))



# ==== RegionRemoval =================
Spectrum_data.RR <- RegionRemoval(Spectrum_data.B, typeofspectra = "serum")

par(default.par) 
plot(Re(Spectrum_data.RR[spectrIndex,]), col = col, xaxt = "n",
     type = "l", ylab = "Intensity", xlab = "ppm", 
     main = "Spectrum after Region Removal")
axis(side = 1, at = at, labels = round(ppmvalues[at],2))



# ==== Normalization =================
Spectrum_data.N <- Normalization(Spectrum_data.RR)

par(mar=c(4,4,1,1), mfrow=c(2,1))
plot(Re(Spectrum_data.RR[spectrIndex,]), col = col1, xaxt = "n",
     type = "l", ylab = "Intensity", xlab = "ppm", 
     main = "Spectrum before Normalization")
axis(side = 1, at = at, labels = round(ppmvalues[at],2))
lines(Re(Spectrum_data.RR[2,]), col = "blue")
lines(Re(Spectrum_data.RR[3,]), col = "green")

plot(Re(Spectrum_data.N[spectrIndex,]), col = col1, xaxt = "n",
     type = "l", ylab = "Intensity", xlab = "ppm", 
     main = "Spectrum after Normalization")
axis(side = 1, at = at, labels = round(ppmvalues[at],2))
lines(Re(Spectrum_data.N[2,]), col = "blue")
lines(Re(Spectrum_data.N[3,]), col = "green")



```
