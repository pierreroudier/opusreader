[![CRAN_Status_Badge](http://www.r-pkg.org/badges/version/opusreader)](https://cran.r-project.org/package=opusreader)
[![Total_Downloads](http://cranlogs.r-pkg.org/badges/grand-total/opusreader)](https://cran.r-project.org/package=opusreader)

# opusreader

Reading OPUS Binary Files in R

## Scope

This package implements a driver for FTIR spectroscopy data collected in the OPUS binary format. 

## Acknowledgement

This package is a very simple fork from the `read_opus_bin_univ` from [the `simplerspec` package authored by Philipp Baumann](https://github.com/philipp-baumann/simplerspec), who deserves all of the credit for this contribution. Here the `simplerspec` functions are modified so they can read `raw` streams directly, and to keep the functions in a lighter package.
 
## Installation

You can use the `remotes` package to install `opusreader`:

`devtools::install_github("pierreroudier/opusreader")`
