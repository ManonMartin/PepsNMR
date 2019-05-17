#' @export DrawPCA
DrawPCA <- function(Signal_data, drawNames = TRUE, main = "PCA score plot", 
                    Class = NULL, axes = c(1, 2), type.pca = c("scores", "loadings"), 
                    loadingstype = c("l", "p"), num.stacked = 4, xlab = "rowname", 
                    createWindow){
  
  
  # Data initialisation and checks ----------------------------------------------
  loadingstype <- match.arg(loadingstype)
  type.pca <- match.arg(type.pca)
  
  checkArg(main, "str", can.be.null = TRUE)
  

  
  # Class
  if (!is.null(Class) & length(Class) != nrow(Signal_data)) {
    stop("the length of class is not equal to the nrow of Signal_data")
  }
  
  
  if (!is.null(Class)) {
    Class_factor <- as.factor(Class)
    nameClass <- deparse(substitute(Class))
  }
  
  
  # axes
  if (!is.vector(axes, mode = "numeric")) {
    stop("axes is not a numeric vector")
  }
  
  if (nrow(Signal_data) < 2) {
    stop("At least 2 spectra are needed for PCA.")
  }
  
  # if (0 %in% Class) {
  #   Class <- Class + 1
  # }
  
  # axes for scores plot
  Xax <- axes[1]
  Yax <- axes[2]
  
  
  # PCA ----------------------------------------------
  
  pca <- stats::prcomp(Re(Signal_data))
  
  # Eigenvalues
  eig <- (pca$sdev)^2
  
  # Variances in percentage
  variance <- eig * 100/sum(eig)
  
  # scores
  scores <- as.data.frame(pca$x)
  
  # loadings
  loadings <- matrix(data = NA, nrow = nrow(pca$rotation), ncol = ncol(pca$rotation))
  for (i in 1:length(eig)) {
    loadings[, i] <- pca$rotation[, i] * pca$sdev[i]
  }
  rownames(loadings) <- colnames(Signal_data)
  colnames(loadings) <- paste0("Loading", c(1:length(eig)))
  loadings <- as.data.frame(loadings)
  
  
  
  # Drawing ----------------------------------------------
  
  plots <- list()
  
  Var <- rowname <- value <- NULL  # only for R CMD check
  
  # SCORES ===============
  if (type.pca == "scores") {
    
    if (!length(axes) == 2) {
      stop("the length of axes is not equal to 2 for scores plot")
    }
    
    if (createWindow)  {
      grDevices::dev.new(noRStudioGD = TRUE)
    }
    Xlim <- c(min(pca$x[, Xax]) * 1.4, max(pca$x[, Xax]) * 1.4)
    Ylim <- c(min(pca$x[, Yax]) * 1.4, max(pca$x[, Yax]) * 1.4)
    
    plots <- ggplot2::ggplot(scores, ggplot2::aes(get(colnames(scores)[Xax]), 
      get(colnames(scores)[Yax]))) + ggplot2::xlim(Xlim) + ggplot2::ylim(Ylim)
    
    if (is.null(Class))  {
      plots <- plots + ggplot2::geom_jitter()
    } else  {
      plots <- plots + ggplot2::geom_point(ggplot2::aes(colour = Class, shape = Class))  +
        ggplot2::scale_shape_discrete(name = nameClass, breaks = unique(Class_factor),
                                      labels = as.character(unique(Class)),
                                      guide=guide_legend(order=1)) +
        ggplot2::scale_colour_discrete(name = nameClass, breaks = unique(Class_factor),
                              labels = as.character(unique(Class)),
                              guide=guide_legend(order=1))
      
    }
    
    plots <- plots + ggplot2::ggtitle(main) + ggplot2::geom_vline(xintercept = 0, 
      size = 0.1) + ggplot2::geom_hline(yintercept = 0, size = 0.1) + 
      ggplot2::theme_bw() + ggplot2::theme(panel.grid.major = ggplot2::element_line(color = "gray60", 
      size = 0.2), panel.grid.minor = ggplot2::element_blank(), panel.background = ggplot2::element_rect(fill = "gray98")) + 
      ggplot2::labs(x = paste0("PC", Xax, " (", round(variance[Xax], 
        2), "%)"), y = paste0("PC", Yax, " (", round(variance[Yax], 
        2), "%)"))
    
    if (drawNames)  {
      if (is.null(Class))  {
        plots <- plots + ggplot2::geom_text(ggplot2::aes(x = scores[, 
          Xax], y = scores[, Yax], label = rownames(Signal_data)), 
          hjust = 0, nudge_x = (Xlim[2]/25), show.legend = FALSE, 
          size = 2)
      } else  {
        plots <- plots + ggplot2::geom_text(ggplot2::aes(x = scores[, 
          Xax], y = scores[, Yax], label = rownames(Signal_data), 
          colour = Class), hjust = 0, nudge_x = (Xlim[2]/25), 
          show.legend = FALSE, size = 2)
      }
    }
    
    print(ggplot2::last_plot())
    
    
  } else {
    # LOADINGS ===============
    
    loadings <- loadings[, axes]
    
    if (is.vector(loadings)) {
      n <- 1
    } else {
      n <- ncol(loadings)
    }
    
    
    i <- 1
    while (i <= n)  {
      if (createWindow)  {
        grDevices::dev.new(noRStudioGD = TRUE)
      }
      
      last <- min(i + num.stacked - 1, n)
      
      if (n == 1)   {
        melted <- reshape2::melt(t(loadings), varnames = c("rowname",  "Var"))
      } else {
        melted <- reshape2::melt(t(loadings[, i:last]), varnames = c("rowname", 
          "Var"))
      }
      
      plot_data <- data.frame(
        facet = paste0("L", i:last," ",
                       paste0("(", round(variance[i:last],2), "%)")),
        x = melted$Var,
        y = melted$value
      )
      
      # to silence CMD check
      # in the future you will be able to use ggplot2::vars(.data$facet)
      # but this is currently not supported (tidyverse/ggplot2#2963)
      facet <- NULL; rm(facet)
      
      plots <- ggplot2::ggplot(data = plot_data, ggplot2::aes(x = .data$x, y = .data$y))
      
      if (loadingstype == "l")  {
        plots <- plots + ggplot2::geom_line()
      } else if (loadingstype == "p")  {
        plots <- plots + ggplot2::geom_point(size = 0.5)
      } else  {
        warning("loadingstype is misspecified")
      }
      
      plots <-  plots +   ggplot2::facet_grid(facet ~  ., scales = "free_y") 
      
      plots <- plots + ggplot2::ggtitle(main) + 
        ggplot2::theme(legend.position = "none") + 
        ggplot2::labs(x = xlab, y = "Loadings") + 
        ggplot2::geom_hline(yintercept = 0, size = 0.5, linetype = "dashed", colour = "gray60") 
      
      if ((melted[1, "Var"] - melted[(dim(melted)[1]), "Var"]) > 0)  {
        plots <- plots + ggplot2::scale_x_reverse()
      }
      
      # Plot finalisation ----------------------------------------------
      
      i <- last + 1
      print(ggplot2::last_plot())
    }
    
  }
  
  
}
