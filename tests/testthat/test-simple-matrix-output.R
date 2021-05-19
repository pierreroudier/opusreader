test_that("simple matrix output works", {

  fn <- opus_file()

  s1 <- opus_read(rep(fn, 5), type = "spec", simplify = TRUE)
  s2 <- opus_read(rep(fn, 5), type = "sc_sample", simplify = TRUE)
  s3 <- opus_read(rep(fn, 5), type = "sc_ref", simplify = TRUE)

  expect_equal(dim(s1$spectra), c(5, 4819))
  expect_equal(dim(s2$spectra), c(5, 4819))
  expect_equal(dim(s3$spectra), c(5, 4819))

})
