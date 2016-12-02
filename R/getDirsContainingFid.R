getDirsContainingFid <- function(path) {
  subdirs <- dir(path, full.names = TRUE)
  if (length(subdirs) > 0) {
    cond <- sapply(subdirs, function(x)  {
      content <- dir(x)
      # subdirs must contain fid, acqu and acqus files
      return("fid" %in% content && "acqu" %in% content && "acqus" %in% content)
    })
    subdirs <- subdirs[cond]
  }
  return(subdirs)
}
