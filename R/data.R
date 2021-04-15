#' @include opusreader.R

#' @name test_spectra.0
#' @title Sample OPUS file containing an MIR soil spectrum
#' @description Sample OPUS binary file. See Details section for more information.
#' @docType data
#' @format A binary OPUS file
#' @details The spectrum contained in this OPUS file was collected
#' on a soil sample from New Zealand using a Brucker FTIR spectrometer.
#'
#' @examples
#' # Access the location of the ASD file using the following command
#' fn <- opus_file()
#' fn
#' # This function is actually just a shorthand for
#' fn <- system.file("extdata", "test_spectra.0", package = "opusreader")
#' fn
NULL

#' @name opus_file
#' @title Get location of a sample OPUS file
#' @description Utility function that retrieves the location of the sample OPUS binary file on disc.
#' @return a character vector storing the location of the sample OPUS file
#' @export
#' @examples
#' fn <- opus_file()
#' fn
opus_file <- function() {
  system.file("extdata", "test_spectra.0", package = "opusreader")
}
