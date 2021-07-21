<!-- NEWS.md is maintained by https://cynkra.github.io/fledge, do not edit -->

# opusreader 0.4.2

- added workaround for corrupted date fields (spits NA in this case)
- Fixed the extra FXV/LXV correction method
- Implemented new spectra types naming scheme (https://github.com/pierreroudier/opusreader/issues/2)
- Changed option name `extract` to `type` (https://github.com/pierreroudier/opusreader/issues/2)
- Added test for simplified output
- implemented correct resampling using wavenumbers based on spectra type (#11)
- Removed occurences of `sapply` with `vapply` (#10)
- added continuous integration using Github Actions


# opusreader 0.4.1

- Fixed the extra FXV/LXV correction method


# opusreader 0.4.0

- Implemented new spectra types naming scheme (https://github.com/pierreroudier/opusreader/issues/2)
- Changed option name `extract` to `type` (https://github.com/pierreroudier/opusreader/issues/2)
- Added test for simplified output
- implemented correct resampling using wavenumbers based on spectra type (#11)
- Removed occurrences of `sapply` with `vapply` (#10)
- added continuous integration using Github Actions


# opusreader 0.3.1.9000

* Added a `NEWS.md` file to track changes to the package.
