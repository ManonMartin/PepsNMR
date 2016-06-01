#' @export DrawPCA
DrawPCA <- function (Signal_data, drawNames=TRUE, createWindow=F, main = "PCA score plot", class = NULL, axes =c(1,2)) {
  checkArg(main, "str", can.be.null=TRUE)
  
  # class
  if(!is.vector(class, mode = "numeric") & !is.null(class)){
    stop("class is not a numeric vector")
  }
  if(is.vector(class, mode = "numeric") & length(class)!=nrow(Signal_data)){
    stop("the length of class is not equal to the nrow of Signal_data")
  }
  
  # axes
  if(!is.vector(axes, mode = "numeric")){
    stop("axes is not a numeric vector")
  }
  
  if(length(axes)!=2){
    stop("the length of axes is not equal to 2")
  }
  
  if(0 %in%class){
    class = class+1
  }
  
  Xax = axes[1]
  Yax = axes[2]
  
  pca <- prcomp(Re(Signal_data))
  
  # Eigenvalues
  eig <- (pca$sdev)^2
  
  # Variances in percentage
  variance <- eig*100/sum(eig)
  
  if (nrow(Signal_data) < 2) {
    stop("At least 2 spectra are needed for PCA.")
  }
  if (createWindow) {
    dev.new(noRStudioGD = TRUE)
  }
  Xlim=c(min(pca$x[,Xax])*1.4, max(pca$x[,Xax])*1.4)
  Ylim=c(min(pca$x[,Yax])*1.4, max(pca$x[,Yax])*1.4)
  
  
  if(is.null(class)) {
    plot(pca$x[,Xax], pca$x[,Yax], xlab="PC1", ylab="PC2")
    if (drawNames) {
      text(pca$x[,Xax], pca$x[,Yax], rownames(Signal_data), pos=c(2,3))
    }
  }else{
    
    
    plot(pca$x[,Xax], pca$x[,Yax], col=class,
         xlab=paste0("PC",Xax," (", round(variance[Xax],2) ,"%)"), xlim=Xlim,
         ylab=paste0("PC",Yax," (", round(variance[Yax],2) ,"%)"), ylim=Ylim,
         main=main)
    
    if (drawNames) {
      text(pca$x[,Xax],pca$x[,Yax],labels=rownames(Signal_data), pos=c(2,3), col=class)
    }
  }
  
  
}
