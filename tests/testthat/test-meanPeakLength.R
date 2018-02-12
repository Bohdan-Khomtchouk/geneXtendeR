context("meanPeakLength and plot")

if (file.exists("peaks.txt")) {
  file.remove("peaks.txt")
}

start <- 500
end <- 1500
by <- 100

test_that("correct message is thrown before peaksinput is run", {
  expect_error(meanPeakLength(rat, start, end), message = "Please run peaksInput() function first!  See ?peaksInput for more information")
  expect_error(meanPeakLengthPlot(rat, start, end, by), message = "Please run peaksInput() function first!  See ?peaksInput for more information")
})

peaksInput(system.file("extdata", "somepeaksfile.txt", package="geneXtendeR"))

test_that("meanPeakLength works properly", {
  expect_equal(meanPeakLength(rat, start, start), 0)
  expect_equal(meanPeakLength(rat, start, end), 3390.304, tolerance = .002, scale = 1)
})

test_that("meanPeakLengthPlot also works properly", {
  expect_silent(meanPeakLengthPlot(rat, start, end, by))
})

test_that("cleanup is successful for meanPeakLength functions", {
  expect_true(file.exists("peaks.txt"))
  expect_silent(file.remove("peaks.txt", sprintf("geneXtender_gtf_%s.bed", seq(start, end, by))))
})