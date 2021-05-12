#' @title Read Bruker OPUS Spectrum Binary Files
#'
#' @description
#' Read file(s) acquired with a Bruker Vertex FTIR Instrument.
#'
#' @param file Character vector with path to file(s).
#' @param extract Character vector of spectra types to extract from OPUS binary file.
#'   Default is `"spc"`, which will extract the final spectra, e.g. expressed in absorbance
#'   (named `AB` in Bruker OPUS program). Possible  additional values for the character vector
#'   supplied to extract are `"spc_nocomp"` (spectrum of the sample without background correction),
#'   `"ScSm"` (single channel spectrum of the sample measurement), `ScRf` (single channel spectrum
#'   of the reference measurement), `"IgSm"` (interferogram of the sample measurement) and `"IgRf"`
#'   (interferogram of the reference measurement).
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
#'     - `metadata`: A data.frame containing metadata from the OPUS file.
#'     - `spc`: If `extract = "spc"`, a matrix of the spectrum of the sample (otherwise set to `NULL`).
#'     - `spc_nocomp`: If `extract = "spc_nocomp"`, a matrix of the spectrum of the sample without
#'       background correction (otherwise set to `NULL`).
#'     - `sc_sm`: If `extract = "ScSm"`, a matrix with the single channel spectrum of the sample
#'       (otherwise set to `NULL`).
#'     - `sc_rf`: If `extract = "ScRf"`, a matrix with the single channel spectrum of the reference
#'       (otherwise set to `NULL`).
#'     - `ig_sm`: If `extract = "IgSm"`, a matrix with the interferogram of the sample
#'       (otherwise set to `NULL`).
#'     - `ig_rf`: If `extract = "IgRf"`, a matrix with the interferogram of the reference
#'       (otherwise set to `NULL`).
#'     - `wavenumbers`:  If `extract = "spc"`, a numeric vector with the wavenumbers of the sample
#'       spectrum (otherwise set to `NULL`).
#'     - `wavenumbers_sc_sm`: If `extract = "ScSm"`, a numeric vector of the wavenumbers of the
#'       single channel spectrum of the sample (otherwise set to `NULL`).
#'     - `wavenumbers_sc_rf`: If `extract = "ScRf"`, a numeric vector with the wavenumbers of the
#'       single channel spectrum of the reference (otherwise set to `NULL`).
#'
#'  - If `simplify = TRUE`, a list of two elements is returned:
#'     - `wavenumbers`: Numeric vector with wavenumbers of the requested spectra.
#'     - `spectra`: Matrix with spectra of requested type (see argument `extract`).
#'
#' @include opus_read_raw.R
#' @importFrom stats approx
#' @export
#'
#' @author Philipp Baumann
#'
opus_read <- function(
  file,
  extract = "spc",
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
        out <- opus_read_raw(rw, extract = extract, atm_comp_minus4offset = atm_comp_minus4offset)

        return(out)
      })
  } else {
    lapply(
      file,
      function(fn) {

        if (!file.exists(fn)) stop(paste0("File '", fn, "' does not exist"))

        # Get raw vector
        rw <- readBin(fn, "raw", 1e5)
        out <- opus_read_raw(rw, extract = extract, atm_comp_minus4offset = atm_comp_minus4offset)

        return(out)
      })
  }

  # If there was only one file to read, we un-nest the list one level
  if (length(file) == 1) {
    res <- res[[1]]
  } else {
    # If a simplified output ( = spectra matrix) was requested
    if (simplify) {

      if (length(extract) > 1) {
        stop("
             Simple output is currently only implemented for one value of the extract option.\n
             A workaround this limitation is to use the `lapply` function, e.g.:\n\n
             lapply(c('spc', 'ScRf'), function(x) read_opus(file, extract = x, simplify = TRUE))
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

          id <- switch(extract,
            spc = "spc",
            spc_nocomp = "spc_nocomp",
            ScSm = "sc_sm",
            ScRf = "sc_rf",
            IgSm = "ig_sm",
            IgRf = "ig_rf"
          )

          # Linear interpolation to get spectra at rounded wavenumber
          s <- approx(x = x$wavenumbers, y = x[[id]], xout = wn_ref, method = "linear")$y
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
