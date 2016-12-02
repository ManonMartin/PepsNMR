#' @export Draw
Draw <- function(Signal_data, type.draw = c("signal", "pca"), output = c("default", 
                "window", "png", "pdf"), dirpath = ".", filename = "%003d", height = 480, 
                width = 640, pdf.onefile = TRUE, ...) {
  
  # Data initialisation and checks ----------------------------------------------
  type.draw <- match.arg(type.draw)
  output <- match.arg(output)
  fullpath <- paste(file.path(dirpath, filename), output, sep = ".")
  createFile <- TRUE
  createWindow <- FALSE
  
  # Drawing --------------------------------------------------------------------
  # output
  switch(output, default = {
    createFile <- FALSE
  }, window = {
    createWindow <- TRUE
    createFile <- FALSE
  }, png = {
    grDevices::png(fullpath, width, height)
  }, pdf = {
    grDevices::pdf(fullpath, width = width/72, height = height/72, 
                  onefile = pdf.onefile)
  }, {
    stop("Unknown output type.")
  })
  
  # Drawing type (signal/spectrum or PCA)
  funs <- list(signal = DrawSignal, pca = DrawPCA)
  if (type.draw %in% names(funs)) {
    fun <- funs[[type.draw]]
  } else {
    stop(paste("Unknown type:", type.draw))
  }
  
  # Plot finalisation ----------------------------------------------
  if (is.vector(Signal_data)) {
    Signal_data <- vec2mat(Signal_data)
  }
  fun(Signal_data, createWindow = createWindow, ...)
  if (createFile) {
    grDevices::dev.off()
  }
}
