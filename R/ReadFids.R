#' @export ReadFids
ReadFids <- function(path, l = 1, subdirs = FALSE) {
  
  # Data initialisation and checks ----------------------------------------------
  begin_info <- beginTreatment("ReadFids")
  checkArg(path, c("str"))
  checkArg(l, c("pos"))
  if (file.exists(path) == FALSE) {
    stop(paste("Invalid path:", path))
  }
  
  
  # Extract the FIDs and their info ----------------------------------------------
  
  if (subdirs == FALSE) {
    fidDirs <- getDirsContainingFid(path)
    n <- length(fidDirs)
    if (n == 0L)  {
      stop(paste("No valid fid in", path))
    }
    fidNames <- sapply(X = fidDirs, FUN = getTitle, l = l, USE.NAMES = F)
    for (i in 1:n)  {
      fidList <- ReadFid(fidDirs[i])
      fid <- fidList[["fid"]]
      info <- fidList[["params"]]
      m <- length(fid)
      if (i == 1)  {
        Fid_data <- matrix(nrow = n, ncol = m, dimnames = list(fidNames, 
          info[["DT"]] * (0:(m - 1))))
        Fid_info <- matrix(nrow = n, ncol = length(info), dimnames = list(fidNames, 
          names(info)))
      }
      Fid_data[i, ] <- fid
      Fid_info[i, ] <- unlist(info)
    }
    
  } else {
    maindirs <- dir(path, full.names = TRUE)
    Fid_data <- numeric()
    Fid_info <- numeric()
    
    for (j in maindirs) {
      fidDirs <- getDirsContainingFid(j)
      n <- length(fidDirs)
      if (n == 0L)  {
        stop(paste("No valid fid in", path))
      }
      fidNames <- sapply(X = fidDirs, FUN = getTitle, l = l, subdirs = subdirs, USE.NAMES = F)
      for (i in 1:n)  {
        fidList <- ReadFid(fidDirs[i])
        fid <- fidList[["fid"]]
        info <- fidList[["params"]]
        m <- length(fid)
        if (i == 1)  {
          Fid_DATA <- matrix(nrow = n, ncol = m, dimnames = list(fidNames, 
          info[["DT"]] * (0:(m - 1))))
          Fid_INFO <- matrix(nrow = n, ncol = length(info), dimnames = list(fidNames, 
          names(info)))
        }
        Fid_DATA[i, ] <- fid
        Fid_INFO[i, ] <- unlist(info)
      }
      Fid_data <- rbind(Fid_data, Fid_DATA)
      Fid_info <- rbind(Fid_info, Fid_INFO)
    }
    
  }
  
  # Check for non-unique IDs ----------------------------------------------
  NonnuniqueIds <- sum(duplicated(row.names(Fid_data)))
  cat("dim Fid_data: ", dim(Fid_data), "\n")
  cat("IDs: ", rownames(Fid_data), "\n")
  cat("non-unique IDs?", NonnuniqueIds, "\n")
  if (NonnuniqueIds > 0) {
    warning("There are duplicated IDs: ", Fid_data[duplicated(Fid_data)])
  }
  
  
  # Return the results ----------------------------------------------
  return(list(Fid_data = endTreatment("ReadFids", begin_info, Fid_data), Fid_info = Fid_info))
  
}
