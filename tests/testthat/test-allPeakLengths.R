context("All peaks length tester")

test_that("all peaks length graphs produces properly", {
  peaksGraph <- allPeakLengths(fpath <- system.file("extdata", "somepeaksfile.txt", package="geneXtendeR"))
  
  expect_equal(peaksGraph$n, 25089)
})

if (file.exists("peaks.txt")) {
  file.remove("peaks.txt")
}