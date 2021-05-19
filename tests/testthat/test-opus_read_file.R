test_that("reading example file works", {

  fn <- opus_file()
  spec <- opus_read(file = fn, progress = FALSE)

  expect_equal(dim(spec$spec), c(1, 4819))
})

test_that("reading different spectral elements works", {

  fn <- opus_file()

  spec <- opus_read(file = fn, type = "spec", progress = FALSE)
  expect_true(!is.null(spec$spec))

  spec <- opus_read(file = fn, type = "sc_sample", progress = FALSE)
  expect_true(!is.null(spec$sc_sample))

  spec <- opus_read(file = fn, type = "sc_ref", progress = FALSE)
  expect_true(!is.null(spec$sc_ref))

  # Support of Bruker's way to name things
  spec <- opus_read(file = fn, type = "ScSm", progress = FALSE)
  expect_true(!is.null(spec$sc_sample))

  spec <- opus_read(file = fn, type = "ScRf", progress = FALSE)
  expect_true(!is.null(spec$sc_ref))

})
