#' @export ReadFids
ReadFids <- function(path, l = 1, subdirs = FALSE, dirs.names = FALSE) {
  
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
    if (dirs.names) {
      separator <- .Platform$file.sep
      path_elem <- strsplit(fidDirs,separator)
      fidNames <- sapply(path_elem, function(x) x[[length(path_elem[[1]])]])
    }else {fidNames <- sapply(X = fidDirs, FUN = getTitle, l = l, subdirs = subdirs,  USE.NAMES = F)}
    
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
    
  } else  {
    maindirs <- dir(path, full.names = TRUE) # subdirectories
    Fid_data <- numeric()
    Fid_info <- numeric()
    
    fidDirs <- c()
    for (j in maindirs) {
      fd <- getDirsContainingFid(j) # recoved FIDs from subdirectories
      n <- length(fd)
      if (n > 0L)  {
        fidDirs <- c(fidDirs, fd)
      } else {warning(paste("No valid fid in",j ))}
    }
     
    if (dirs.names==TRUE) {
      if (length(fidDirs)!= length(dir(path))) { # at least one subdir contains more than 1 FID
        separator <- .Platform$file.sep
        path_elem <- strsplit(fidDirs,separator)
        fidNames <- sapply(path_elem, function(x) paste(x[[length(path_elem[[1]])-1]],
                                                        x[[length(path_elem[[1]])]], sep = "_"))
      }else {fidNames <- dir(path)}
      
    } else {fidNames <- sapply(X = fidDirs, FUN = getTitle, l = l, subdirs = subdirs, USE.NAMES = F)}
    
    for (i in 1:length(fidNames))  {
      fidList <- ReadFid(fidDirs[i])
      fid <- fidList[["fid"]]
      info <- fidList[["params"]]
      m <- length(fid)
      if (i == 1)  {
        Fid_data <- matrix(nrow = length(fidNames), ncol = m, dimnames = list(fidNames, 
                                                                              info[["DT"]] * (0:(m - 1))))
        Fid_info <- matrix(nrow = length(fidNames), ncol = length(info), dimnames = list(fidNames, 
                                                                                         names(info)))
      }
      Fid_data[i, ] <- fid
      Fid_info[i, ] <- unlist(info)
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
