context(".geneXtender checker") #hidden method

fpath <- system.file("extdata", "somepeaksfile.txt", package="geneXtendeR")
peaksInput(fpath)

test_that("test dimensions and make sure it works", {
  expect_equal(dim(.geneXtender(2500, rat, TRUE)), c(32402, 5))
})

if (file.exists("peaks.txt")) {
  file.remove("peaks.txt")
}