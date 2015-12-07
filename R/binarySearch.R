binarySearch <- function (a, target, lower=TRUE) {
  # search the index i in a such that a[i] == target
  # if it doesn't exists and  lower, it searches the closer a[i] such that a[i] < target
  #                          !lower, it seraches the closer a[i] such that a[i] > target
  # a should be monotone but can be increasing or decreasing

  # if a is increasing
  # INVARIANT: a[amin] < target < a[amax]
  N <- length(a)
  if ((a[N] - target) * (a[N] - a[1]) <= 0) {
    return(N)
  }
  if ((a[1] - target) * (a[N] - a[1]) >= 0) {
    return(1)
  }
  amin <- 1
  amax <- N
  while (amin + 1 < amax) {
    amid <- floor((amin + amax) / 2)
    if ((a[amid] - target) * (a[amax] - a[amid]) < 0) {
      amin <- amid
    } else if ((a[amid] - target) * (a[amax] - a[amid]) > 0) {
      amax <- amid
    } else {
      # a[amid] == a[amax] or a[amid] == target
      # In both cases, a[amid] == target
      return(amid)
    }
  }
  if (xor(lower, a[amin] > a[amax])) {
    # (lower && a[amin] < a[amax]) || (!lower && a[min] > a[max])
    # If increasing and we want the lower, the take amin
    # If decreasing and we want the bigger, we take amin too
    return(amin)
  } else {
    return(amax)
  }
}