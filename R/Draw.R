Draw <- function (Signal_data, type = c("signal","pca"),
                  output = c("default","window","png","pdf"), dirpath = ".",
                  filename = "%003d", height = 480, width = 640, pdf.onefile = TRUE, ...) {
  type = match.arg(type)
  output = match.arg(output)
  fullpath = paste(file.path(dirpath, filename), output, sep = ".")
  createFile <- TRUE
  createWindow <- FALSE
  switch(output,
         "default" = {
           createFile <- FALSE
         },
         "window" = {
           createWindow <- TRUE
           createFile <- FALSE
         },
         "png" = {
           png(fullpath, width, height)
         },
         "pdf" = {
           pdf(fullpath, width = width/72, height = height/72, onefile = pdf.onefile)
         },
         {
            stop('Unknown output type.')
         }
  )
  funs <- list(signal=DrawSignal, pca=DrawPCA)
  if (type %in% names(funs)) {
    fun <- funs[[type]]
  } else {
    stop(paste("Unknown type:", type))
  }

  if (is.vector(Signal_data)) {
    Signal_data = vec2mat(Signal_data)
  }
  fun(Signal_data, createWindow=createWindow, ...)
  if (createFile) {
    dev.off()
  }
}
