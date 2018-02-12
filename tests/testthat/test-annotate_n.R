context("annotate_n checker")

peaksInput(system.file("extdata", "somepeaksfile.txt", package="geneXtendeR"))
n <- 5
extension <- 1234

test_that("annotate_n is the correct length", {
  expect_silent(annotated <- annotate_n(rat, extension, n))
  annotated <- annotate_n(rat, extension, n)
  
  expect_equal(dim(annotated)[1], dim(samplepeaksinput)[1] * n)
  
  expect_equal(dim(annotated)[2], 11)
  
  expect_equal(names(annotated), c("Peak-Num", "Chromosome", "Peak-Start", "Peak-End", "Gene-Start", "Gene-End", "Gene-ID", "Gene-Name", "rank", "Minimum-Distance-to-Gene", "seqid"))
})


test_that("cleanup is successful for annotate_n", {
  expect_silent(file.remove("peaks.txt"))
  expect_true(file.exists(sprintf("annotated_%s_%s.txt", extension, n)))
  expect_silent(file.remove(sprintf("annotated_%s_%s.txt", extension, n)))
  expect_true(file.exists(sprintf("geneXtender_gtf_%s.bed", extension)))
  expect_silent(file.remove(sprintf("geneXtender_gtf_%s.bed", extension)))
})