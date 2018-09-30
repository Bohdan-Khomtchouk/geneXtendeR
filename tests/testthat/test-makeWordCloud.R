context("makeWordCloud tester")
library(org.Rn.eg.db)

if (file.exists("peaks.txt")) {
  file.remove("peaks.txt")
}

start <- 500
end <- 1500
GOspecies <- org.Rn.eg.db

test_that("correct message is thrown before peaksinput is run", {
  expect_message(makeWordCloud(rat, start, end, BP, GOspecies), message = "Please run peaksInput() function first!  See ?peaksInput for more information")
})

peaksInput(system.file("extdata", "somepeaksfile.txt", package="geneXtendeR"))

test_that("error is thrown when not selecting correct category", {
  expect_error(makeWordCloud(rat, start, end, AP, GOspecies), "Not a valid GO category.  Must be either BP, CC, or MF.")
})

test_that("cleanup is successful for makeNetwork", {
  expect_silent(file.remove("peaks.txt"))
  expect_silent(file.remove(sprintf("geneXtender_gtf_%s.bed", c(start, end))))
})