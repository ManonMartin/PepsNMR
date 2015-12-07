DrawPCA <- function (Signal_data, drawNames=TRUE, createWindow=F) {
  pca <- prcomp(Re(Signal_data))
  if (nrow(Signal_data) < 2) {
    stop("At least 2 spectra are needed for PCA.")
  }
  if (createWindow) {
    dev.new(noRStudioGD = TRUE)
  }
  plot(pca$x[,1], pca$x[,2], xlab="PC1", ylab="PC2")
  if (drawNames) {
    text(pca$x[,1], pca$x[,2], rownames(Signal_data), pos=3)
  }
}