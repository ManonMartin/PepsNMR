#' @export ReadFids
ReadFids <- function(path, l=1) {
  begin_info <- beginTreatment("ReadFids")
  checkArg(path, c("str"))
  checkArg(l, c("pos"))
  fidDirs <- getDirsContainingFid(path)
  n <- length(fidDirs)
  if (n == 0L) {
    stop(paste("No valid fid in", path))
  }
  fidNames <- sapply(X = fidDirs, FUN = getTitle, l = l, USE.NAMES=F)
  for (i in 1:n) {
    fidList <- ReadFid(fidDirs[i])
    fid <- fidList[["fid"]]
    info <- fidList[["params"]]
    m <- length(fid)
    if (i == 1) {
      Fid_data <- matrix(nrow=n, ncol=m, dimnames=list(fidNames, info[["DT"]] * (0:(m-1))))
      Fid_info <- matrix(nrow=n, ncol=length(info), dimnames=list(fidNames,names(info)))
    }
    Fid_data[i,] = fid
    Fid_info[i,] = unlist(info)
  }
  return(list(Fid_data=endTreatment("ReadFids", begin_info, Fid_data), Fid_info=Fid_info))
}
