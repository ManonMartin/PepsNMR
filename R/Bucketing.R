#' @export Bucketing
Bucketing <- function (Spectrum_data, m = 500) {
  left_part_trapz <- function (x1, x2, xmid, y1, y2) {
    # Integral from x1 to xmid of the trapezium defined by (x1,y1) and (x1,y2)
    #        /| y2
    #       / |
    #      /  |          # : Integral
    #     /#  |
    # y1 /##  |
    #    |##  |
    #    |##  |
    #    x1 | x2
    #      xmid
    if (x1 == x2) {
      # Integral is 0 but should be of the form of y1 and y2
      return(y1 * 0)
    } else {
      ymid <- y1 + (y2 - y1) * (xmid - x1) / (x2 - x1)
      half_interval = abs(xmid - x1) / 2
      return((y1 + ymid) * half_interval)
    }
  }
  begin_info <- beginTreatment("Bucketing", Spectrum_data)
  Spectrum_data <- begin_info[["Signal_data"]]
  checkArg(m, "int", "pos")
  ppm <- as.numeric(colnames(Spectrum_data))
  n <- nrow(Spectrum_data)
  old_m <- ncol(Spectrum_data)
  if (n == 0) {
    stop("Empty ppm scale")
  }
  decreasing = (ppm[old_m] <= ppm[1])
  buckets <- seq(ppm[1], ppm[old_m], length.out=m+1)
  centers <- (buckets[1:m] + buckets[2:(m+1)]) / 2
  bucketed <- matrix(0, nrow=n, ncol=m, dimnames=list(rownames(Spectrum_data), centers))
  for (i in 1:m) {
    # ppm[from] is at the left of buckets[i].
    # if ppm it is decreasing, that means we want it higher e.g. ppm[from] >= buckets[i]
    from = binarySearch(ppm, buckets[i], !decreasing)
    to = binarySearch(ppm, buckets[i+1], decreasing)

    # Now that ppm[from] is on the left and ppm[to] is on the right,
    # we can take the trapezium over this interval and remove
    #  * the integral from ppm[from] to buckets[i] and
    #  * the integral from buckets[i+1¯ to ppm[to].
    if (from < to && buckets[i] != buckets[i+1]) { # else the integral is just 0
      # Trapezium from ppm[from] to ppm[to]] (we will need to remove some at the extremities)
      half_intervals_length <- abs(ppm[(from+1):to] - ppm[from:(to-1)]) / 2
      # The contribution of the points to half the area of the trapezium on their right
      trapz_right <- base::rowSums(sweep(Spectrum_data[, from:(to-1), drop=F], MARGIN=2, half_intervals_length, `*`))
      # The contribution of the points to half the area of the trapezium on their left
      trapz_left  <- base::rowSums(sweep(Spectrum_data[, (from+1):to, drop=F], MARGIN=2, half_intervals_length, `*`))
      # Integral of the part of the integral of the leftmost trapezium that is after buckets[i]
      left_part   <- left_part_trapz(ppm[from], ppm[from+1], buckets[i]  , Spectrum_data[, from], Spectrum_data[, from+1])
      right_part  <- left_part_trapz(ppm[to]  , ppm[to-1]  , buckets[i+1], Spectrum_data[, to], Spectrum_data[, to-1])
      bucketed[, i] = (trapz_right + trapz_left - left_part - right_part) / (abs(buckets[i+1] - buckets[i]))
    }
  }
  return(endTreatment("Bucketing", begin_info, bucketed))
}
