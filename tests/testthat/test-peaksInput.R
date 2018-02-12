context("peaks-input checker")

test_that("peaksInput inputs and null", {
  fpath <- system.file("extdata", "somepeaksfile.txt", package="geneXtendeR")
  peaksInput(fpath)
  
  expect_error(peaksInput(''))
  expect_true(file.exists("peaks.txt"))
  file.remove("peaks.txt")
})

if (file.exists("peaks.txt")) {
  file.remove("peaks.txt")
}
