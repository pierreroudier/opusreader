#' @title Read a Bruker OPUS spectrum binary file
#'
#' @description
#' Read single binary file acquired with an
#' Bruker Vertex FTIR Instrument
#'
#' @param file Character vector with path to file
#' @param extract Character vector of spectra types to extract from OPUS binary
#' file. Default is \code{"spc"}, which will extract the final spectra, e.g.
#' expressed in absorbance (named \code{AB} in Bruker OPUS programs). Possible
#' additional values for the character vector supplied to extract are \code{"spc_nocomp"} (spectrum of the sample without background correction),
#' \code{"ScSm"} (single channel spectrum of the sample measurement), \
#' code{"ScRf"} (single channel spectrum of the reference measurment),
#' \code{"IgSm"} (interferogram of the sample measurment) and \code{"IgRf"}
#' (interferogram of the reference measurement).
#' @param simplify Logical (defaults \code{FALSE}): if set to \code{TRUE}, returns a much smaller list. The first object of that list (\code{wavenumbers}) is the wavenumbers of the first file read. The second object (\code{spectra}) is a matrix of the corresponding spectra. Especially useful when passing more than one file to the \code{file} option, for example to read a suite of spectral file directly into a matrix.
#' @param wns_digits Integer that specifies the number of decimal places used to round
#' the wavenumbers (values of x-variables) if \code{simplify = TRUE}.
#' @param progress Logical (defaults \code{TRUE}) whether a message is
#' printed when an OPUS binary file is parsed into an R list entry.
#' @param atm_comp_minus4offset Logical whether spectra after atmospheric
#' compensation are read with an offset of \code{-4} bites from Bruker OPUS
#' files. Default is \code{FALSE}.
#'
#' @return The output of \code{opus_read} depends on the value of the \code{simplify} option used in the function call.
#'
#'  - if \code{simplify = FALSE} (default), \code{opus_read} returns a list of 10 elements:
#'     - \code{metadata}: a \code{data.frame} containing metadata from the OPUS file
#'     - \code{spc} If \code{"spc"} was requested in the \code{extract} option, a matrix of the spectrum of the sample (otherwise set to \code{NULL}).
#'     - \code{spc_nocomp} If \code{"spc_nocomp"} was requested in the \code{extract} option, a matrix of the spectrum of the sample without background correction (otherwise set to \code{NULL}).
#'     - \code{sc_sm} If \code{"ScSm"} was requested in the \code{extract} option, a matrix of the single channel spectrum of the sample (otherwise set to \code{NULL}).
#'     - \code{sc_rf} If \code{"ScRf"} was requested in the \code{extract} option, a matrix of the single channel spectrum of the reference (otherwise set to \code{NULL}).
#'     - \code{ig_sm} If \code{"IgSm"} was requested in the \code{extract} option, a matrix of the interferogram of the sample (otherwise set to \code{NULL}).
#'     - \code{ig_rf}  If \code{"IgRf"} was requested in the \code{extract} option, a matrix of the interferogram of the reference (otherwise set to \code{NULL}).
#'     - \code{wavenumbers} If \code{"spc"} was requested in the \code{extract} option, a numeric vector of the wavenumbers of the spectrum of the sample (otherwise set to \code{NULL}).
#'     - \code{wavenumbers_sc_sm} If \code{"ScSm"} was requested in the \code{extract} option, a numeric vector of the wavenumbers of the single channel spectrum of the sample (otherwise set to \code{NULL}).
#'     - \code{wavenumbers_sc_rf} If \code{"ScRf"} was requested in the \code{extract} option, a numeric vector of the wavenumbers of the single channel spectrum of the reference (otherwise set to \code{NULL}).
#'
#'  - if \code{simplify = TRUE}:
#'     - \code{wavenumbers} numeric, wavenumbers of the spectrum requested.
#'     - \code{spectra} matrix, requested spectrum.
#'
#' @usage opus_read(file, extract = "spc", simplify = FALSE, wns_digits = 1L, progress = TRUE,
#' atm_comp_minus4offset = FALSE)
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
        rw <- readBin(fn, "raw", 10e9)
        out <- opus_read_raw(rw, extract = extract, atm_comp_minus4offset = atm_comp_minus4offset)

        return(out)
      })
  } else {
    lapply(
      file,
      function(fn) {

        if (!file.exists(fn)) stop(paste0("File '", fn, "' does not exist"))

        # Get raw vector
        rw <- readBin(fn, "raw", 10e9)
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

          if (extract == "spc") id <- "spc"
          else if (extract == "ScSm") id <- "sc_sm"
          else if (extract == "ScRf") id <- "sc_rf"
          else if (extract == "IgSm") id <- "ig_sm"
          else if (extract == "IgRf") id <- "ig_rf"

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
