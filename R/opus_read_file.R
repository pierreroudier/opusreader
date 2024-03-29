#' @title Read Bruker OPUS Spectrum Binary Files
#'
#' @description
#' Read file(s) acquired with a Bruker Vertex FTIR Instrument.
#'
#' @param file Character vector with path to file(s).
#' @param type Character vector of spectra types to extract from OPUS binary
#' file. Default is `"spec"`, which will extract the final spectra, e.g.
#' expressed in absorbance (named `AB` in Bruker OPUS programs). Possible
#' additional values for the character vector supplied to `type` are `"spec_no_atm_comp"` (spectrum of the sample without compensation for atmospheric gases, water vapor and/or carbon dioxide),
#' `"sc_sample"` (single channel spectrum of the sample measurement), `"sc_ref"` (single channel spectrum of the reference measurement),
#' `"ig_sample"` (interferogram of the sample measurement) and `"ig_ref"`
#' (interferogram of the reference measurement).
#'
#' @param simplify Logical (defaults to `FALSE`): if set to `TRUE`, returns a flattened list.
#'   The first element of that list (`wavenumbers`) is the wavenumbers of the first file read.
#'   The second element (`spectra`) is a matrix of the corresponding spectra. Especially useful when
#'   passing more than one file to the `file` option, for example to read a suite of spectral file.
#'   directly into a matrix.
#'   If the files do not have the same number of wavebands, the wavebands of the very first file passed to `file` is used as a reference for linear interpolation.
#' @param wns_digits Integer that specifies the number of decimal places used to round
#'   the wavenumbers (values of x-variables) if `simplify = TRUE`.
#' @param progress Logical (defaults to `TRUE`) whether a message is printed when an OPUS binary file
#'   is parsed into an R list entry.
#' @param atm_comp_minus4offset Logical whether spectra after atmospheric compensation are read with
#'   an offset of -4 bytes from Bruker OPUS files. Default is `FALSE`.
#'
#' @details The type of spectra returned by the function when using `type = "spec"` depends on the setting of the Bruker instrument: typically, it can be either absorbance or reflectance.
#'
#' The type of spectra to extract from the file can also use Bruker's OPUS software naming conventions, as follows:
#'
#' - `ScSm` corresponds to `sc_sample`
#' - `ScRf` corresponds to `sc_ref`
#' - `IgSm` corresponds to `ig_sample`
#' - `IgRf` corresponds to `ig_ref`
#'
#' @return The output of `opus_read()` depends on the value of the `simplify` option used in the function call.
#'
#'  - If `simplify = FALSE` (default), `opus_read()` returns a list of 10 elements:
#'     - `metadata`: a `data.frame` containing metadata from the OPUS file
#'     - `spec` If `"spec"` was requested in the `type` option, a matrix of the spectrum of the sample (otherwise set to `NULL`).
#'     - `spec_no_atm_comp` If `"spec_no_atm_comp"` was requested in the `type` option, a matrix of the spectrum of the sample without background compensation (otherwise set to `NULL`).
#'     - `sc_sample` If `"sc_sample"` was requested in the `type` option, a matrix of the single channel spectrum of the sample (otherwise set to `NULL`).
#'     - `sc_ref` If `"sc_ref"` was requested in the `type` option, a matrix of the single channel spectrum of the reference (otherwise set to `NULL`).
#'     - `ig_sample` If `"ig_sample"` was requested in the `type` option, a matrix of the interferogram of the sample (otherwise set to `NULL`).
#'     - `ig_ref`  If `"ig_ref"` was requested in the `type` option, a matrix of the interferogram of the reference (otherwise set to `NULL`).
#'     - `wavenumbers` If `"spec"` or `"spec_no_atm_comp"` was requested in the `type` option, a numeric vector of the wavenumbers of the spectrum of the sample (otherwise set to `NULL`).
#'     - `wavenumbers_sc_sample` If `"sc_sample"` was requested in the `type` option, a numeric vector of the wavenumbers of the single channel spectrum of the sample (otherwise set to `NULL`).
#'     - `wavenumbers_sc_ref` If `"sc_ref"` was requested in the `type` option, a numeric vector of the wavenumbers of the single channel spectrum of the reference (otherwise set to `NULL`).

#'
#'  - If `simplify = TRUE`, a list of two elements is returned:
#'     - `wavenumbers`: Numeric vector with wavenumbers of the requested spectra.
#'     - `spectra`: Matrix with spectra of requested type (see argument `type`).
#'
#' @include opus_read_raw.R
#' @importFrom stats approx
#' @export
#'
#' @author Philipp Baumann
#'
opus_read <- function(
  file,
  type = "spec",
  simplify = FALSE,
  wns_digits = 1L,
  progress = TRUE,
  atm_comp_minus4offset = FALSE
) {

  res <- if (requireNamespace("pbapply", quietly = TRUE) & progress) {
    pbapply::pblapply(
      file,
      function(fn) {

        if (!file.exists(fn)) stop(paste0("File '", fn, "' does not exist"))

        # Get raw vector
        rw <- readBin(fn, "raw", n = file.size(fn))
        out <- opus_read_raw(rw, type = type, atm_comp_minus4offset = atm_comp_minus4offset)

        return(out)
      })
  } else {
    lapply(
      file,
      function(fn) {

        if (!file.exists(fn)) stop(paste0("File '", fn, "' does not exist"))

        # Get raw vector
        rw <- readBin(fn, "raw", n = file.size(fn))
        out <- opus_read_raw(rw, type = type, atm_comp_minus4offset = atm_comp_minus4offset)

        return(out)
      })
  }

  # If there was only one file to read, we un-nest the list one level
  if (length(file) == 1) {
    res <- res[[1]]
  } else {
    # If a simplified output ( = spectra matrix) was requested
    if (simplify) {

      if (length(type) > 1) {
        stop("
             Simple output is currently only implemented for one value of the `type` option.\n
             A workaround this limitation is to use the `lapply` function, e.g.:\n\n
             lapply(c('spec', 'sc_ref'), function(x) read_opus(file, type = x, simplify = TRUE))
             ", call. = FALSE)
      }

      # Fetch wavenumbers
      wns <- lapply(
        res,
        function(x) round(x$wavenumbers, digits = wns_digits)
      )

      # Arbitrarily take the first rounded WN as the reference one
      wn_ref <- wns[[1]]

      # Check the wavenumbers have all the same length
      if (length(unique(sapply(wns, length))) > 1) {
        # stop("Spectra can't be combined since they don't all have the same number of wavenumbers.", call. = FALSE)
        warning("Spectra don't all have the same number of wavenumbers. Interpolation will be used to combine them in a matrix.", call. = FALSE)
      }

      specs <- lapply(
        res,
        function(x) {

          id <- switch(type,
            spec = "spec",
            spec_no_atm_comp = "spec_no_atm_comp",
            sc_sample = "sc_sample",
            sc_ref = "sc_ref",
            ig_sample = "ig_sample",
            ig_ref = "ig_ref"
          )

          # Grab correct wavenumbers for interpolation
          wn <- switch(type,
            spec = x$wavenumbers,
            spec_no_atm_comp = x$wavenumbers,
            sc_sample = x$wavenumbers_sc_sample,
            sc_ref = x$wavenumbers_sc_ref,
            ig_sample = x$wavenumbers_sc_sample,
            ig_ref = x$wavenumbers_sc_ref
          )

          # Linear interpolation to get spectra at rounded wavenumber
          s <- approx(
            x = wn,
            y = x[[id]],
            xout = wn_ref,
            method = "linear"
          )$y
          names(s) <- as.character(wn_ref)

          return(s)
        })

      res <- list(
        wavenumbers = wns[[1]],
        spectra = do.call(rbind, specs)
      )

    }
  }

  return(res)
}
