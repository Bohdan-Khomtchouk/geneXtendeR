context("Makenetwork tester")
library(org.Rn.eg.db)

if (file.exists("peaks.txt")) {
  file.remove("peaks.txt")
}

start <- 500
end <- 1500
GOspecies <- org.Rn.eg.db

test_that("correct message is thrown before peaksinput is run", {
  expect_message(makeNetwork(rat, start, end, BP, GOspecies), message = "Please run peaksInput() function first!  See ?peaksInput for more information")
})

peaksInput(system.file("extdata", "somepeaksfile.txt", package="geneXtendeR"))

test_that("error is thrown when not selecting correct category", {
  expect_error(makeNetwork(rat, start, end, AP, GOspecies), "Not a valid GO category.  Must be either BP, CC, or MF.")
})

test_that("makeNetwork runs properly", {
  expect_message(makeNetwork(rat, start, end, BP, GOspecies))
  expect_message(makeNetwork(rat, start, end, CC, GOspecies))
  expect_message(makeNetwork(rat, start, end, MF, GOspecies))
})

test_that("cleanup is successful for makeNetwork", {
  expect_silent(file.remove("peaks.txt"))
  expect_silent(file.remove(sprintf("geneXtender_gtf_%s.bed", c(start, end))))
})