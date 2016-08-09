#' @export DrawPCA
DrawPCA <- function(Signal_data, drawNames=TRUE, main = "PCA score plot", Class = NULL, axes =c(1,2),
                     type =c("scores", "loadings"), loadingstype=c("l", "p"), num.stacked=4, xlab="rowname", createWindow) {
  
  loadingstype=match.arg(loadingstype)
  type = match.arg(type)
  
  checkArg(main, "str", can.be.null=TRUE)
  
  # Class
  if(!is.vector(Class, mode = "any") & !is.null(Class)){
    stop("Class is not a numeric vector")
  }
  if(is.vector(Class, mode = "numeric") & length(Class)!=nrow(Signal_data)){
    stop("the length of Class is not equal to the nrow of Signal_data")
  }
  
  # axes
  if(!is.vector(axes, mode = "numeric")){
    stop("axes is not a numeric vector")
  }
  
  if(length(axes)<2){
    stop("the length of axes is not equal or superior to 2")
  }
  
  if (nrow(Signal_data) < 2) {
    stop("At least 2 spectra are needed for PCA.")
  }
  
  if(0 %in%Class){
    Class = Class+1
  }
  
  Xax = axes[1]
  Yax = axes[2]
  
  pca <- stats::prcomp(Re(Signal_data))
  
  # Eigenvalues
  eig <- (pca$sdev)^2
  
  # Variances in percentage
  variance <- eig*100/sum(eig)
  
  # scores
  scores = as.data.frame(pca$x)
  
  # loadings
  loadings = matrix(data = NA, nrow = nrow(pca$rotation), ncol = ncol(pca$rotation))
  for (i in 1:length(eig)) {
    loadings[,i] = pca$rotation[,i]*pca$sdev[i]
  }
  rownames(loadings) = colnames(Signal_data)
  colnames(loadings) = paste0("Loading", c(1:length(eig)))
  loadings = as.data.frame(loadings)
  
  plots <- list()
  
  if (type == "scores") {
    
 
  if (createWindow) {
    grDevices::dev.new(noRStudioGD = TRUE)
  }
  Xlim=c(min(pca$x[,Xax])*1.4, max(pca$x[,Xax])*1.4)
  Ylim=c(min(pca$x[,Yax])*1.4, max(pca$x[,Yax])*1.4)
    
    plots <- ggplot2::ggplot(scores, ggplot2::aes(get(colnames(scores)[Xax]),get(colnames(scores)[Yax]))) +
      ggplot2::xlim(Xlim) +
      ggplot2::ylim(Ylim) 
      
      if(is.null(Class)) {
        plots <- plots + ggplot2::geom_jitter()
      } else {plots <- plots +  ggplot2::geom_jitter(ggplot2::aes(colour = Class, shape = Class))}
      
      plots <- plots + ggplot2::ggtitle(main) +
      ggplot2::geom_vline(xintercept = 0, size = 0.1) +
      ggplot2::geom_hline(yintercept = 0, size = 0.1) +
      ggplot2::theme_bw() +
      ggplot2::theme(panel.grid.major = ggplot2::element_line(color = "gray60", size = 0.2), panel.grid.minor = ggplot2::element_blank(), 
                     panel.background = ggplot2::element_rect(fill = "gray98")) +
      ggplot2::labs(x=paste0("PC",Xax," (", round(variance[Xax],2) ,"%)"), y=paste0("PC",Yax," (", round(variance[Yax],2) ,"%)")) 
      
      if (drawNames) {  
        if(is.null(Class)) {
          plots = plots + ggplot2::geom_text(ggplot2::aes(x = scores[,Xax], y = scores[,Yax], label = rownames(Signal_data)),  
                                              hjust = 0, nudge_x = (Xlim[2]/25),  show.legend = FALSE, size = 2)
        } else {plots = plots + ggplot2::geom_text(ggplot2::aes(x = scores[,Xax], y = scores[,Yax], label = rownames(Signal_data), colour = Class, shape = Class),  
                                                    hjust = 0, nudge_x = (Xlim[2]/25), show.legend = F, size = 2)}
      }

  print(ggplot2::last_plot())
  
  } else {
    loadings = loadings[,axes]
    n = ncol(loadings)
    
    i = 1
    while (i <= n)
    {
      if (createWindow) {
        grDevices::dev.new(noRStudioGD = TRUE) 
      }
     
      last = min(i + num.stacked-1, n)
      
      melted <- reshape2::melt(t(loadings[, i:last]), varnames=c("rowname", "Var"))
      
      plots <- ggplot2::ggplot(data = melted, ggplot2::aes(x = Var, y = value)) 
      if (loadingstype == "l") {
        plots = plots + ggplot2::geom_line() 
        
      } else if (loadingstype == "p") {
        
        plots = plots + ggplot2::geom_point(size = 0.5) 
      } else {warning("loadingstype is misspecified")}
      
      plots = plots + ggplot2::ggtitle(main) +
        ggplot2::facet_grid(rowname ~ ., scales = "free_y") +
        ggplot2::theme(legend.position="none") +
        ggplot2::labs(x=xlab, y = "Loadings") +
        ggplot2::geom_hline(yintercept = 0, size = 0.5, linetype = "dashed", colour = "gray60") +
        ggplot2::annotate("text", x = -Inf, y = Inf, label = paste0("(",round(variance[i:last],2), "%)"), vjust=1, hjust=1)
      
      if ((melted[1,"Var"] - melted[(dim(melted)[1]),"Var"])>0) {
        plots =  plots + ggplot2::scale_x_reverse() 
      }
      # 
      #         require("gridExtra")
      i = last + 1
      print(ggplot2::last_plot())
    }
    
  }
  
  
}
