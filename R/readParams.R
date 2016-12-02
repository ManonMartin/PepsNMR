# Read parameter values for Fid_info in the ReadFids function

readParams <- function(file, paramsName) {
  
  isDigit <- function(c) {
    return(suppressWarnings(!is.na(as.numeric(c))))
  }
  lines <- readLines(file)
  params <- sapply(paramsName, function(x) NULL)
  
  for (paramName in paramsName)  {
    # Find the line with the parameter I add a '$' '=' in the pattern so that for
    # example 'TD0' is not found where I look for 'TD' and LOCSW and WBSW when I look
    # for 'SW'
    pattern <- paste("\\$", paramName, "=", sep = "")
    occurences <- grep(pattern, lines)
    if (length(occurences) == 0L)  {
      stop(paste(file, "has no field", pattern))
    }
    if (length(occurences) > 1L) {
      warning(paste(file, "has more that one field", pattern, " I take the first one"))
    }
    line <- lines[occurences[1]]
    
    # Cut beginning and end of the line '##$TD= 65536' -> '65536'
    first <- 1
    while (first <= nchar(line) & !isDigit(substr(line, first, first))) {
      first <- first + 1
    }
    last <- nchar(line)
    while (last > 0 & !isDigit(substr(line, last, last)))  {
      last <- last - 1
    }
    params[paramName] <- as.numeric(substr(line, first, last))
  }
  return(params)
}
