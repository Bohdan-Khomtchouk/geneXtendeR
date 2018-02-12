context("annotate()")

if (file.exists("peaks.txt")) {
  file.remove("peaks.txt")
}
extension <- 1234

test_that("peaks file needs to be created first", {
  expect_error(annotate(rat, extension), message = "Please run peaksInput() function first!  See ?peaksInput for more information")
})

peaksInput(system.file("extdata", "somepeaksfile.txt", package="geneXtendeR"))

test_that("annotated file exists and is correct", {
  annotated <- annotate(rat, extension)
  
  expect_true(file.exists(sprintf("peaks_annotated_%s.txt", extension)))
  expect_equal(names(annotated), c("Chromosome", "Peak-Start", "Peak-End", "Gene-Start", "Gene-End", "Gene-ID", "Gene-Name", "Distance-of-Gene-to-Nearest-Peak"))
  expect_equal(dim(annotated)[1], dim(samplepeaksinput)[1])
})

test_that("cleanup for this method works", {
  expect_silent(file.remove("peaks.txt"))
  expect_silent(file.remove(sprintf("peaks_annotated_%s.txt", extension)))
  expect_true(file.exists(sprintf("geneXtender_gtf_%s.bed", extension)))
  expect_silent(file.remove(sprintf("geneXtender_gtf_%s.bed", extension)))
})