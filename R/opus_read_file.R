#' @title Read Bruker OPUS Spectrum Binary Files
#'
#' @description
#' Read file(s) acquired with a Bruker Vertex FTIR Instrument.
#'
#' @param file Character vector with path to file(s).
#' @param type Character vector of spectra types to extract from OPUS binary
#' file. Default is \code{"spec"}, which will extract the final spectra, e.g.
#' expressed in absorbance (named \code{AB} in Bruker OPUS programs). Possible
#' additional values for the character vector supplied to \code{type} are \code{"spec_no_bc"} (spectrum of the sample without background correction),
#' \code{"sc_sample"} (single channel spectrum of the sample measurement), \
#' code{"sc_ref"} (single channel spectrum of the reference measurement),
#' \code{"ig_sample"} (interferogram of the sample measurement) and \code{"ig_ref"}
#' (interferogram of the reference measurement).
#'
#' @param simplify Logical (defaults to `FALSE`): if set to `TRUE`, returns a flattened list.
#'   The first element of that list (`wavenumbers`) is the wavenumbers of the first file read.
#'   The second element (`spectra`) is a matrix of the corresponding spectra. Especially useful when
#'   passing more than one file to the `file` option, for example to read a suite of spectral file
#'   directly into a matrix.
#' @param wns_digits Integer that specifies the number of decimal places used to round
#'   the wavenumbers (values of x-variables) if `simplify = TRUE`.
#' @param progress Logical (defaults to `TRUE`) whether a message is printed when an OPUS binary file
#'   is parsed into an R list entry.
#' @param atm_comp_minus4offset Logical whether spectra after atmospheric compensation are read with
#'   an offset of `-4` bytes from Bruker OPUS files. Default is `FALSE`.
#'
#' @return The output of `opus_read()` depends on the value of the `simplify` option used in the function call.
#'
#'  - If `simplify = FALSE` (default), `opus_read()` returns a list of 10 elements:
#'     - \code{metadata}: a \code{data.frame} containing metadata from the OPUS file
#'     - \code{spec} If \code{"spec"} was requested in the \code{type} option, a matrix of the spectrum of the sample (otherwise set to \code{NULL}).
#'     - \code{spec_no_bc} If \code{"spec_no_bc"} was requested in the \code{type} option, a matrix of the spectrum of the sample without background correction (otherwise set to \code{NULL}).
#'     - \code{sc_sample} If \code{"sc_sample"} was requested in the \code{type} option, a matrix of the single channel spectrum of the sample (otherwise set to \code{NULL}).
#'     - \code{sc_ref} If \code{"sc_ref"} was requested in the \code{type} option, a matrix of the single channel spectrum of the reference (otherwise set to \code{NULL}).
#'     - \code{ig_sample} If \code{"ig_sample"} was requested in the \code{type} option, a matrix of the interferogram of the sample (otherwise set to \code{NULL}).
#'     - \code{ig_ref}  If \code{"ig_ref"} was requested in the \code{type} option, a matrix of the interferogram of the reference (otherwise set to \code{NULL}).
#'     - \code{wavenumbers} If \code{"spec"} or \code{"spec_no_bc"} was requested in the \code{type} option, a numeric vector of the wavenumbers of the spectrum of the sample (otherwise set to \code{NULL}).
#'     - \code{wavenumbers_sc_sample} If \code{"sc_sample"} was requested in the \code{type} option, a numeric vector of the wavenumbers of the single channel spectrum of the sample (otherwise set to \code{NULL}).
#'     - \code{wavenumbers_sc_ref} If \code{"sc_ref"} was requested in the \code{type} option, a numeric vector of the wavenumbers of the single channel spectrum of the reference (otherwise set to \code{NULL}).

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
        rw <- readBin(fn, "raw", 1e5)
        out <- opus_read_raw(rw, type = type, atm_comp_minus4offset = atm_comp_minus4offset)

        return(out)
      })
  } else {
    lapply(
      file,
      function(fn) {

        if (!file.exists(fn)) stop(paste0("File '", fn, "' does not exist"))

        # Get raw vector
        rw <- readBin(fn, "raw", 1e5)
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
             ")
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
        stop("Spectra can't be combined since they don't all have the same number of wavenumbers.")
      }

      specs <- lapply(
        res,
        function(x) {

          id <- switch(type,
            spec = "spec",
            spec_no_bc = "spec_no_bc",
            sc_sample = "sc_sample",
            sc_ref = "sc_ref",
            ig_sample = "ig_sample",
            ig_ref = "ig_ref"
          )

          # Grab correct wavenumbers for interpolation
          wn <- switch(type,
            spec = x$wavenumbers,
            spec_no_bc = x$wavenumbers,
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
