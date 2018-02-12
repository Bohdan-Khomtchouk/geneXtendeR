context("diffGo tests")

if (file.exists("peaks.txt")) {
  file.remove("peaks.txt")
}

start <- 500
end <- 1500
GOspecies <- org.Rn.eg.db

test_that("correct message is thrown before peaksinput is run", {
  expect_message(diffGO(rat, start, end, BP, GOspecies), message = "Please run peaksInput() function first!  See ?peaksInput for more information")
})

peaksInput(system.file("extdata", "somepeaksfile.txt", package="geneXtendeR"))

test_that("error is thrown when not selecting correct category", {
  expect_message(diffGO(rat, start, end, AP, GOspecies), "Not a valid GO category.  Must be either BP, CC, or MF.")
})

test_that("function runs properly for each category", {
  expect_message(p <- diffGO(rat, start, end, BP, GOspecies))
  expect_message(p <- diffGO(rat, start, end, CC, GOspecies))
  expect_message(p <- diffGO(rat, start, end, MF, GOspecies))
})

test_that("function output is correct for BP", {
  p <- diffGO(rat, start, end, BP, GOspecies)
  expect_equal(dim(p), c(683, 3))
  expect_length(unique(p$`gene$SYMBOL`), 75)
})


test_that("cleanup is successful for diffGO", {
  expect_true(file.exists("peaks.txt"))
  expect_silent(file.remove(sprintf("geneXtender_gtf_%s.bed", c(start, end)), "peaks.txt"))
})
