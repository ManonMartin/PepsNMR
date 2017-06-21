# Get the name of the signal from the title file or fromt the name of the subdirectory

getTitle <- function(path, l, subdirs) {
  title <- NULL
  title_file <- file.path(file.path(file.path(path, "pdata"), "1"), "title")
  if (file.exists(title_file)) {
    lines <- readLines(title_file, warn = FALSE)
    if (length(lines) >= 1)  {
      first_line <- gsub("^\\s+|\\s+$", "", lines[l])
      if (nchar(first_line) >= 1)  {
        title <- first_line
      } else {
        warning(paste("The", l ,"line of the title file is blank for directory ", 
          path, "and the (sub)dirs names are used instead"))
      }
    } else {
      warning(paste("The title file is empty for directory ", path, "and the (sub)dirs names are used instead"))
    }
  } else {
    warning(paste("Title file doesn't exists for directory ", path, "\n the (sub)dirs names are  used instead"))
  }
  if (is.null(title)) {
    if(subdirs) {
      separator <- .Platform$file.sep
      path_elem <- strsplit(path,separator)[[1]]
      title <- paste(path_elem[length(path_elem)-1], path_elem[length(path_elem)], sep = "_")
    } else{title <- basename(path)} 
  }
  return(title)
}
