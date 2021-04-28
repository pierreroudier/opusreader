test_that("reading example file works", {

  fn <- opus_file()
  spec <- opus_read(file = fn, progress = FALSE)

  expect_equal(dim(spec$spc), c(1, 4819))
})

test_that("reading different spectral elements works", {

  fn <- opus_file()

  spec <- opus_read(file = fn, extract = "spc", progress = FALSE)
  expect_true(!is.null(spec$spc))

  spec <- opus_read(file = fn, extract = "ScSm", progress = FALSE)
  expect_true(!is.null(spec$sc_sm))

  spec <- opus_read(file = fn, extract = "ScRf", progress = FALSE)
  expect_true(!is.null(spec$sc_rf))

})
