context("Hotspotplot")

allpeaks <- system.file("extdata", "totalpeaksfile.txt", package="geneXtendeR")
sigpeaks <- system.file("extdata", "significantpeaksfile.txt", package="geneXtendeR")
start <- 0
end <- 5000
by <- 200

test_that("hotspotPlot works", {
  expect_silent(hotspotPlot(allpeaks, sigpeaks, rat, start, end, by))
})

test_that("cleanup is successful for hotspotPlot", {
  expect_silent(file.remove("peaks.txt"))
  expect_silent(file.remove(sprintf("geneXtender_gtf_%s.bed", seq(start, end, by))))
})