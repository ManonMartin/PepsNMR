#' @export DrawSignal
DrawSignal <-
function (Signal_data,
                        subtype=c("together", "separate", "stacked", "diffmean", "diffmedian", "diffwith"),
                        ReImModArg=c(T, F, F, F),
                        vertical=T,
                        xlab="rowname",
                        main.names=NULL,
                        nticks=42,
                        row=1,
                        num.stacked=4,
                        createWindow=F) {
  subtype  <- match.arg(subtype)

  scale    <- colnames(Signal_data)
  n        <- nrow(Signal_data)
  m        <- ncol(Signal_data)
  vticks   <- round(seq(1, m, length=nticks))
  vlabels  <- scale[vticks]
  vlabels  <- round(as.numeric(vlabels),2)
  num.plot <- sum(ReImModArg)
  if (num.plot <= 0) {
    stop("Nothing selected in ReImModArg.")
  } else if (num.plot <= 2) {
    if (vertical) {
      nrow <- num.plot
      ncol <- 1
    } else {
      nrow <- 1
      ncol <- num.plot
    }
  } else {
    nrow <- 2
    ncol <- 2
  }

  if (is.null(main.names)) {
    main.names <- rownames(Signal_data)
    if (is.null(main.names)) {
      main.names <- 1:n
    }
  } else {
    if (!is.vector(main.names)) {
      stop("main.names is not a vector")
    }
    if (length(main.names) != n) {
      stop(paste("main.names has length", length(main.names), "and there are", n, "FIDs."))
    }
  }

  elements <- list()
  if (ReImModArg[1]) {
    elements[["Re"]]  <- Re(Signal_data)
  }
  if (ReImModArg[2]) {
    elements[["Im"]]  <- Im(Signal_data)
  }
  if (ReImModArg[3]) {
    elements[["Mod"]] <- Mod(Signal_data)
  }
  if (ReImModArg[4]) {
    elements[["Arg"]] <- Arg(Signal_data)
  }
  if (subtype == "separate" | subtype == "stacked") {
    i = 1
    while (i <= n)
    {
      if (createWindow) {
        dev.new(noRStudioGD = TRUE) 
      }
      if (subtype == "separate") {
        # The other uses gridExtra to do that
        par(mfrow=c(nrow,ncol))
      }
      plots <- list()
      if (subtype == "separate") {
        last = i
      } else {
        last = min(i + num.stacked-1, n)
      }
      for (name in names(elements)) {
        if (subtype == "separate") {
          plot(elements[[name]][i,],type="l",main=main.names[i],xaxt="n",ylab=name,xlab=xlab)
          axis(side=1,at=vticks,labels=vlabels,cex.axis=0.6,las=2)
        } else {
#           require("ggplot2")
#           require("reshape2")
          melted <- reshape2::melt(elements[[name]][i:last,], varnames=c("rowname", "Var"))
          plots[[name]] <- ggplot2::ggplot(melted, aes("Var", "value")) +
                  ggplot2::geom_line() +
                  ggplot2::facet_grid(rowname ~ ., scales = "free_y") +
                  ggplot2::theme(legend.position="none") +
                  ggplot2::labs(x=xlab, y=name)
        }
      }
      if (subtype == "stacked") {
#         require("gridExtra")
        do.call(grid.arrange, c(plots, list(nrow=nrow, ncol=ncol)))
      }
      i = last + 1
    }
  } else {
    rainbow_colors <- rainbow(n)
    if (createWindow) {
      dev.new(noRStudioGD = TRUE) 
    }
    par(mfrow=c(nrow,ncol))

    # Loop for Re, Im, Mod and Arg
    for (name in names(elements)) {
      # Get this part of the signal
      element <- elements[[name]]
      # Express the signal according to a reference if asked by `subtype'
      if (subtype == "diffmean") {
        element <- sweep(element, MARGIN=2, colMeans(element), `-`)
      } else if (subtype == "diffmedian") {
#         require("matrixStats")
        element <- sweep(element, MARGIN=2, colMedians(element), `-`)
      } else if (subtype == "diffwith") {
        element <- sweep(element, MARGIN=2, element[row, ], `-`)
        if (row == 1 & n > 1) {
          # Since we use plot on the first row
          # and lines on the following, the y scale is calculated at the first row
          # so if the first row is all 0, it causes problems
          tmp <- element[1, ]
          element[1, ] <- element[2, ]
          element[2, ] <- tmp
        }
      }
      for (i in 1:n)
      {
        if (i == 1) {
          # If it is our first plot, we need to set the axes
          plot(element[i,],col=rainbow_colors[i],type="l",main=main.names[i],xaxt="n",ylab=name,xlab=xlab)
          axis(side=1,at=vticks,labels=vlabels,cex.axis=0.6,las=2)
        } else {
          lines(element[i,],col=rainbow_colors[i])
        }
      }
      legend(x="topright",legend=main.names,lty=1,lwd=2, col=rainbow_colors)
    }
  }
}
