context("linePlot")

if (file.exists("peaks.txt")) {
  file.remove("peaks.txt")
}

start <- 500
end <- 2000
by <- 150

test_that("correct message is thrown before peaksinput is run", {
  expect_message(linePlot(rat, start, end, by), message = "Please run peaksInput() function first!  See ?peaksInput for more information")
})

peaksInput(system.file("extdata", "somepeaksfile.txt", package="geneXtendeR"))

test_that("lineplot runs correctly", {
  expect_silent(linePlot(rat, start, end, by))
})

test_that("cleanup is successful for linePlot", {
  expect_silent(file.remove("peaks.txt"))
  expect_silent(file.remove(sprintf("geneXtender_gtf_%s.bed", seq(start, end, by))))
})