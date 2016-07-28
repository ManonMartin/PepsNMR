getTitle <- function(path, l) {
  title <- NULL
  title_file = file.path(file.path(file.path(path, "pdata"), "1"), "title")
  if (file.exists(title_file)) {
    lines <- readLines(title_file, warn = FALSE)
    if (length(lines) >= 1) {
      first_line = gsub("^\\s+|\\s+$", "", lines[l])
      if (nchar(first_line) >= 1) {
        title <- first_line
      } else {
        warning(paste("The first line of the title file is blank for directory ", path))
      }
    } else {
      warning(paste("The title file is empty for directory ", path))
    }
  } else {
    warning(paste("Title file doesn't exists for directory ", path))
  }
  if (is.null(title)) {
    title = basename(path)
  }
  return (title)
}