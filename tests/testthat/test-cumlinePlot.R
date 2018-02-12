context("CumlinePlot tests")

if (file.exists("peaks.txt")) {
  file.remove("peaks.txt")
}

start <- 500
end <- 1500
by <- 100

test_that("correct message is thrown before peaksinput is run", {
  expect_message(cumlinePlot(rat, start, end, by), message = "Please run peaksInput() function first!  See ?peaksInput for more information")
})

peaksInput(system.file("extdata", "somepeaksfile.txt", package="geneXtendeR"))

test_that("cumlinePlot works correctly", {
  expect_silent(cumlinePlot(rat, start, end, by))
})

test_that("cleanup is successful for cumlinePlot", {
  expect_true(file.exists("peaks.txt"))
  expect_silent(file.remove(sprintf("geneXtender_gtf_%s.bed", seq(start, end, by)), "peaks.txt"))
})