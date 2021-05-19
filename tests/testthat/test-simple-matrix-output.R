test_that("simple matrix output works", {

  fn <- opus_file()

  s1 <- opus_read(rep(fn, 5), type = "spc", simplify = TRUE)
  s2 <- opus_read(rep(fn, 5), type = "ScSm", simplify = TRUE)
  s3 <- opus_read(rep(fn, 5), type = "ScRf", simplify = TRUE)

  expect_equal(dim(s1$spectra), c(5, 4819))
  expect_equal(dim(s2$spectra), c(5, 4819))
  expect_equal(dim(s3$spectra), c(5, 4819))

})

