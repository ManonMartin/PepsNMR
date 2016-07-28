require('testthat')

full_test <- function(input, answer, n, old_m, m, rownames, input_colnames, answer_colnames, fliplr) {
  input_Spectrum  <- matrix(input,  nrow=n, ncol=old_m, byrow=T, dimnames=list(rownames, input_colnames))
  answer_Spectrum <- matrix(answer, nrow=n, ncol=m,     byrow=T, dimnames=list(rownames, answer_colnames))
  if (fliplr) {
    input_Spectrum  <- input_Spectrum[, old_m:1, drop=F]
    answer_Spectrum <- answer_Spectrum[, m:1, drop=F]
  }
  output_Spectrum <- Bubble::Bucketing(input_Spectrum, m)
    expect_equal(output_Spectrum, answer_Spectrum)
}

double_full_test <- function(name, input, answer, n, old_m, m, rownames, input_colnames, answer_colnames) {
  test_that(paste("Increasing", name, " "), {
    full_test(input, answer, n, old_m, m, rownames, input_colnames, answer_colnames, F)
  })
  test_that(paste("Decreasing", name, " "), {
    full_test(input, answer, n, old_m, m, rownames, input_colnames, answer_colnames, T)
  })
}

double_full_test("simple as 1,2,3", 1:6, c(2,5), 2, 3, 1, c("a", "b"), c("1","2","3"), c("2"))

double_full_test("simple non equidistant PPM", 1:6, c(2.125,5.125), 2, 3, 1, c("a", "b"), c("0","3","8"), c("4"))

# new ppm             |           |
# borders       |           |           |
# integrals          1.75         4
# ppm           0   1   2   3   4   5   6
# values        3   1     4        -1   5
# interpolated              3
# 7.5 = 1*(3+1)/2 + 1.5*(1+4)/2 + 0.5*(4+3)/2
#   4 = 2*(3-1)/2 + 1*(-1+5)/2
# We divide them by 3 = (3-0) = (6-3)
double_full_test("one harder signal", c(3,1,4,-1,5), c(7.5/3, 4/3), 1, 5, 2, c("pi"), c("0","1","2.5","5","6"), c("1.5","4.5"))