#' @title Read a Bruker OPUS spectrum binary file
#'
#' @description
#' Read single binary file acquired with an
#' Bruker Vertex FTIR Instrument
#'
#' @param file_path Character vector with path to file
#' @param extract Character vector of spectra types to extract from OPUS binary
#' file. Default is \code{c("spc")}, which will extract the final spectra, e.g.
#' expressed in absorbance (named \code{AB} in Bruker OPUS programs). Possible
#' additional values for the character vector supplied to extract are
#' \code{"ScSm"} (single channel spectrum of the sample measurement), \
#' code{"ScRf"} (single channel spectrum of the reference measurment),
#' \code{"IgSm"} (interferogram of the sample measurment) and \code{"IgRf"}
#' (interferogram of the reference measurement).
#' @param print_progress Logical (default \code{TRUE}) whether a message is
#' printed when an OPUS binary file is parsed into an R list entry.
#' @param atm_comp_minus4offset Logical whether spectra after atmospheric
#' compensation are read with an offset of \code{-4} bites from Bruker OPUS
#' files. Default is \code{FALSE}.
#' @usage read_opus_bin_univ(file_path, extract = c("spc"),
#' print_progress = TRUE, atm_comp_minus4offset = FALSE)
# Importing functions `%do%` and foreach::`%dopar%` does not work, see
# http://stackoverflow.com/questions/30216613/how-to-use-dopar-when-only-import-foreach-in-description-of-a-package
# Got the following error:
# "Error : object '`%do%`' is not exported by 'namespace:foreach'"
#' @include opus_read_raw.R
#' @importFrom hexView readRaw
#' @export
#'
#' @author Philipp Baumann
#'
opus_read_file <- function(
  file_path,
  extract = c("spc"),
  print_progress = TRUE,
  atm_comp_minus4offset = FALSE
) {

  if (!file.exists(file_path)) {
    stop(paste0("File does not exist"))
  }

  try({

    # file_path <- "data/soilspec_background/yamsys_bg_gold/BF_lo_15_soil_cal.0"
    # Read entire content of file as bytes

    pa <- hexView::readRaw(file_path, offset = 0,
                           nbytes = file.info(file_path)$size, human = "char",
                           size = 1, endian = "little")

    # Get raw vector
    pr <- pa$fileRaw

    out <- opus_read_raw(pr)

    return(out)
  }) # closes try() function

}
