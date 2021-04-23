#' @title Read a Bruker OPUS spectrum binary file
#'
#' @description
#' Read single binary file acquired with an
#' Bruker Vertex FTIR Instrument
#'
#' @param file_path Character vector with path to file
#' @param extract Character vector of spectra types to extract from OPUS binary
#' file. Default is \code{"spc"}, which will extract the final spectra, e.g.
#' expressed in absorbance (named \code{AB} in Bruker OPUS programs). Possible
#' additional values for the character vector supplied to extract are
#' \code{"ScSm"} (single channel spectrum of the sample measurement), \
#' code{"ScRf"} (single channel spectrum of the reference measurement),
#' \code{"IgSm"} (interferogram of the sample measurement) and \code{"IgRf"}
#' (interferogram of the reference measurement).
#' @param print_progress Logical (default \code{TRUE}) whether a message is
#' printed when an OPUS binary file is parsed into an R list entry.
#' @param atm_comp_minus4offset Logical whether spectra after atmospheric
#' compensation are read with an offset of \code{-4} bites from Bruker OPUS
#' files. Default is \code{FALSE}.
#' @usage opus_read(file_path, extract = "spc",
#' print_progress = TRUE, atm_comp_minus4offset = FALSE)
#' @include opus_read_raw.R
#' @export
#'
#' @author Philipp Baumann
#'
opus_read <- function(
  file_path,
  extract = "spc",
  print_progress = TRUE,
  atm_comp_minus4offset = FALSE
) {

  if (!file.exists(file_path)) {
    stop(paste0("File does not exist"))
  }

  # Get raw vector
  rw <- readBin(file_path, "raw", 10e9)
  out <- opus_read_raw(rw, extract = extract, atm_comp_minus4offset = atm_comp_minus4offset)

  return(out)
}
