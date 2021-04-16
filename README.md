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

`remotes::install_github("pierreroudier/opusreader")`

## Basic usage

The package includes a sample OPUS spectra for demonstration purposes. It slocation can be determined using the function `opus_file()`:

```
> opus_file()
[1] "/path/to/opusreader/inst/extdata/test_spectra.0"
```

Two functions are available to actually read the data:

1. `opus_read` is the workhorse of the package, and allows to read binary files:

```
# File name of the OPUS file to read
fn <- opus_file()

# Read OPUS file
s <- opus_read(fn)

# Plot the spectra
plot(s$wavenumbers, s$spc, type = 'l')
```

2. The `opus_read_raw` function is lower level, and allows to read a `raw` representation of the OPUS file. It is useful if the OPUS file is "streamed" over eg a web service:

```
# Store OPUS file content as raw binary data
fn <- opus_file()
rw <- readBin(fn, "raw", 10e9)

# Read raw binary data directly
s <- opus_read_raw(rw)

# Plot the spectra
plot(s$wavenumbers, s$spc, type = 'l')
```

