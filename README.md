# Tull coude spectrograph pipeline

Reduction and analysis pipeline for the Tull coude spectrograph on the 2.7-m Harlan J. Smith Telescope at McDonald Observatory

[![Powered by Astropy](http://img.shields.io/badge/powered%20by-AstroPy-orange.svg?style=flat)](http://www.astropy.org/)

## Description

This is a end-to-end reduction and analysis pipeline for the Tull coude spectrograph at McDonald Observatory. It starts with raw FITS files, performs basic CCD image processing, extracts 1D spectra, and wavelength calibrates. It also has analysis modules to continuum fit science spectra and measure radial velocities using broadening functions.

This pipeline is based on an original Python 2 reduction pipeline written primarily by Daniel Krolikowski, with initial work by Aaron Rizzuto.
The original pipeline can be found at: https://github.com/dkrolikowski/coudereduction

This new version of the pipeline was written to make it more user friendly and modular, to ease future development or changes. I decided to start fresh rather than make modifications to the old Python 2 pipeline because of the significant changes to the structure and framework of the code.

The pipeline is designed to be modular so that it could be adapted for other instruments or set-ups in the future. Certain modules in this code could be used univerisally and independent of the data source. However, this pipeline is still written with the Tull coude spectrograph in mind, so it has plenty of settings and code that are specific to this instrument.