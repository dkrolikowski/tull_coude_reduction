Generating CCD Calibration Files and Image Processing
=====================================================

CCD Calibration Files
---------------------

- median combine bias and flat field
- generate the bad pixel mask

In the config file, this step's options are set in the ``calibration`` section of the YAML file. 

:doc:`ccd_calibrations`

..
	.. autofunction:: ccd_calibrations.build_calibrations
		:no-index:

CCD Image Processing
--------------------

Steps

- subtract the bias
- flat field
- optional: cosmic ray subtraction