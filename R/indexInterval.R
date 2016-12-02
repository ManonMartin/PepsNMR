# returns the discrete borders of the interval for a numeric vector a

indexInterval <- function (a, from, to, inclusive=TRUE) {
  # If inclusive and from <= to, we need to take the lower
  # If not inclusive and from > to, we need to take the lower too
  lowerFrom <- xor(inclusive, from > to)
  fromIndex <- binarySearch(a, from, lowerFrom)
  toIndex <- binarySearch(a, to, !lowerFrom)
  return(fromIndex:toIndex)
}