# Tull coude spectrograph pipeline

Reduction and analysis pipeline for the Tull coude spectrograph on the 2.7-m Harlan J. Smith Telescope at McDonald Observatory

[![Powered by Astropy](http://img.shields.io/badge/powered%20by-AstroPy-orange.svg?style=flat)](http://www.astropy.org/)

## Description

This is a end-to-end reduction and analysis pipeline for the Tull coude spectrograph at McDonald Observatory. 

For reduction, the pipeline stars with raw FITS files and then:
- Generates CCD calibration files (bias, flat field, bad pixel mask) and perform basic image processing
- Uses the flat field to trace the echelle orders
- Extracts 1D spectra for arc lamps and science observations
- Fits a wavelength solution using arc lamp spectra and wavelength calibrates the science observations

It then proceeds to use analysis modules to:
- Fit a spline continuum to science observations (in this case assuming absorption spectra, like for stars)
- Compute broadening functions to measure stellar radial velocities, in most cases to sub-km/s precision (in the best cases to tens of m/s)

This pipeline is based on an original Python 2 reduction pipeline written primarily by Daniel Krolikowski, with initial work by Aaron Rizzuto.
The original pipeline can be found at: https://github.com/dkrolikowski/coudereduction.

This new version of the pipeline was written to make it more user friendly and modular, to ease future development or changes. I decided to start fresh rather than make modifications to the old Python 2 pipeline because of the significant changes to the structure and framework of the code.

The pipeline is designed to be modular so that it could be adapted for other instruments or set-ups in the future. Certain modules in this code could be used univerisally and independent of the data source. However, this pipeline is still written with the Tull coude spectrograph in mind, so it has plenty of settings and code that are specific to this instrument.

## Documentation

## Installation

## Acknowledgements

This pipeline was primarily developed by [Daniel Krolikowski](https://dkrolikowski.github.io) (Steward Observatory, UT Austin). The first version of the pipeline was developed with input from Dr. Aaron Rizzuto (then UT Austin Postdoc).

I thank Benjamin Tofflemire, Catherine Manea, Adam Kraus, Keith Hawkins, Bill Cochran, and Chris Sneden for conversations, advice, and data sources that enabled the development of this pipeline. Particular thanks to Benjamin Tofflemire as the lead developer of the saphires package used for RV measurement and Shubham Kanodia for writing barycorrpy.

This code depends on:
- [astropy](https://www.astropy.org)
- [astroscrappy](https://astroscrappy.readthedocs.io/en/latest/)
- [barycorrpy](https://github.com/shbhuk/barycorrpy)
- [matplotlib](https://matplotlib.org)
- [numpy](https://numpy.org)
- [pandas](https://pandas.pydata.org)
- [saphires](https://saphires.readthedocs.io/en/latest/index.html)
- [scipy](https://scipy.org)

