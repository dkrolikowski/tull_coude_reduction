Pipeline Overview and Structure
===============================

Overview of the pipeline steps
------------------------------

The pipeline follows the general steps for reducing data from an echelle spectrograph.
These general steps are:

1. Generate a header manifest
2. Gather flats and biases and generate calibration files by combining them. Also BPM
3. Image process: the rest of the images bias subtract, flat field, and cosmic subtract if it says to.
4. Trace the orders! Using the flat field
5. Extract 1D spectra along the echelle order traces. Different methods for different types of spectra.
6. Wavelength calibration: fit a wavelength solution and then apply to the rest of the data.

That gives you the final spectra. We then have some more "analysis" steps:

1. Continuum fit stellar spectra.
2. Calculate radial velocities using broadening functions.

How is the pipeline actually run? Discuss structure: there is a module file for each broad step. Within that module is a function meant to act as the main script to run that step.
Then there is a main reduction run file that wraps all these modules together and takes in a config file that defines the run.

Briefly discuss the config file: it will be touched upon in each page for the modules and then in a separate page where everything is explicitly defined.
There is a config file included in the package that is meant to be a template.

Pipeline structure: how to run it
---------------------------------

I will put instructions on installing the package here once I get to that.

.. _target_to_config_description:

The reduction configuration file
--------------------------------