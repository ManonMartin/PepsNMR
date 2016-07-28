require('testthat')

full_test <- function(input, answer, n, m, rownames, colnames) {
  input_RawSpect  <- matrix(input,  nrow=n, ncol=m, byrow=T, dimnames=list(rownames, colnames))
  answer_RawSpect <- matrix(answer, nrow=n, ncol=m,     byrow=T, dimnames=list(rownames, colnames))
  output_RawSpect <- Bubble::NegativeValuesZeroing(input_RawSpect)
  expect_equal(output_RawSpect, answer_RawSpect)
}

test_that("1 signal", {
  full_test(c(-1,1,-2),
            c( 0,1, 0),
            1, 3, c("marvin"), c("a","x","ax"))
})

test_that("2 signals", {
  full_test(c(-1,1,3,-2),
            c( 0,1,3, 0),
            2, 2, c("mostly","harmless"), c("x","y"))
})

test_that("4 signals", {
  full_test(c(1,2,9,-7,8,-10,.1,-.3,-1,9,0,-9,1,-2,1,3,-9,0,0,-1),
            c(1,2,9, 0,8,  0,.1,  0, 0,9,0, 0,1, 0,1,3, 0,0,0, 0),
            4, 5, c("a","b","c","d"), c("-2", "-1.5", "0", "1", "0"))
})