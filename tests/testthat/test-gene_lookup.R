context("checking function gene_lookup")

if (file.exists("peaks.txt")) {
  file.remove("peaks.txt")
}

extension <- 1234
cutoff <- 0
n <- 5
gene_name <- "Vom2r3"

test_that("correct message is thrown before peaksinput is run", {
  expect_message(gene_lookup(rat, gene_name), message = "Please run peaksInput() function first!  See ?peaksInput for more information")
})

peaksInput(system.file("extdata", "somepeaksfile.txt", package="geneXtendeR"))

test_that("user error handling", {
  expect_message(gene_lookup(rat, gene_name, n = 0, extension, cutoff), message = "Please select a cutoff >= 0 or a value of n >= 1")
})

test_that("gene_lookup works correctly", {
  expect_silent(p <- gene_lookup(rat, gene_name, n, extension))
  p <- gene_lookup(rat, gene_name, n, extension)
  expect_equal(dim(p), c(10, 7))
})

test_that("cutoff works properly", {
  expect_message(gene_lookup(rat, gene_name, n, extension, cutoff = -1), message = "Please select a cutoff >= 0 or a value of n >= 1")
  expect_silent(p <- gene_lookup(rat, gene_name, n, extension, cutoff))
  p <- gene_lookup(rat, gene_name, n, extension, cutoff)
  expect_equal(dim(p), c(0, 7))
})

test_that("File order is chosen properly", {
  .geneXtender(1235, rat, FALSE)
  .geneXtender(0, rat, FALSE)
  p <- gene_lookup(rat, gene_name, n)
  q <- gene_lookup(rat, gene_name, n, extension)
  expect_true(any(p != q, TRUE))
  expect_silent(file.remove(sprintf("geneXtender_gtf_%s.bed", c(1235, 0))))
})

test_that("cleanup is successful for gene_lookup", {
  expect_true(file.exists("peaks.txt"))
  expect_silent(file.remove(sprintf("geneXtender_gtf_%s.bed", extension), "peaks.txt"))
})