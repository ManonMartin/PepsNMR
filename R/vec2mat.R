# vector -> matrix
vec2mat <- function(vec) {
  
  return(matrix(vec, nrow = 1, dimnames = list(c(1), names(vec)))) 
  
  }