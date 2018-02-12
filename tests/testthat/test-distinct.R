context("Distinct tests and peakLengthBoxplot")

if (file.exists("peaks.txt")) {
  file.remove("peaks.txt")
}

start <- 100
end <- 5000

test_that("correct message is thrown before peaksinput is run", {
  expect_error(distinct(rat, start, end), message = "Please run peaksInput() function first!  See ?peaksInput for more information")
  expect_error(peakLegthBoxPlot(rat, start, end), message = "Please run peaksInput() function first!  See ?peaksInput for more information")
})

peaksInput(system.file("extdata", "somepeaksfile.txt", package="geneXtendeR"))

test_that("distinct works properly", {
  expect_silent(distinct(rat, start, end))
  p <- distinct(rat, start, end)
  expect_equal(dim(p), c(1034, 9))
  expect_length(unique(p$Distance), 1)
  expect_silent(peakLengthBoxplot(rat, start, end))
})

test_that("cleanup is successful for cumlinePlot", {
  expect_true(file.exists("peaks.txt"))
  expect_silent(file.remove(sprintf("geneXtender_gtf_%s.bed", c(start, end)), "peaks.txt"))
})