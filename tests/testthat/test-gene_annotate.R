context("gene annotation testing")

if (file.exists("peaks.txt")) {
  file.remove("peaks.txt")
}
extension <- 1234

test_that("correct message is thrown before peaksinput is run", {
  expect_error(gene_annotate(rat, extension), message = "Please run peaksInput() function first!  See ?peaksInput for more information")
})

peaksInput(system.file("extdata", "somepeaksfile.txt", package="geneXtendeR"))

test_that("gene_annotate works properly", {
  expect_silent(gene_annotate(rat, extension))
  p <- gene_annotate(rat, extension)
  expect_equal(dim(p), c(11815, 9))
  expect_equal(sum(p$`Number-of-Peaks-Associated-with-Gene`), dim(data.table::fread("peaks.txt"))[1])
})

test_that("cleanup is successful for gene_annotate", {
  expect_true(file.exists("peaks.txt"))
  expect_silent(file.remove(sprintf("geneXtender_gtf_%s.bed", extension), "peaks.txt"))
})