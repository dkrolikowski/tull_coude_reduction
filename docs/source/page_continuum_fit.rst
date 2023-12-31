Stellar spectrum continuum fit
==============================

After the wavelength calibration step, the "reduction" process of the pipeline is finished. We now have fully extracted and calibrated echelle spectra to do science with!

The first process in the "analysis" portion of the pipeline is to fit the continuum (including remaining blaze) of stellar observations. This step kind of straddles the line between a reduction and analysis process, but we include it in the analyis portion because it is first required for other science measurements, like radial velocities or equivalent widths.

Continuum fitting is performed with the :py:meth:`continuum_fit <modules.continuum_fit>` module. Options are defined in the ``continuum_fit`` section of the main *config* YAML file, which is described in full :ref:`here <target_to_config_description>`.

Fitting the continuum with a spline
-----------------------------------

As written, the pipeline will fit the continuum of all on-sky science observations. For each spectral order, the continuum is fit using the ``scipy.interpolate`` implementation of a B-spline. The fitting is done iteratively, with flux outliers (primarily absorption lines) being excluded from the fit to define only the continuum level.

The *config* file defines the parameters setting the knot spacing of the B-spline. These options include the number of chunks to break the spectrum into and a minimum knot spacing in wavelength (to use if the number of chunks would require a smaller knot spacing).

The iterative fitting is done with the :py:func:`continuum_fit_with_spline <modules.continuum_fit.continuum_fit_with_spline>` function. The *config* file sets the lower and upper sigma rejection levels. For our default, the lower sigma rejection level is lower to preferentially exclude absorption lines. We set the maximum number of iterations to 10.

The continuum spectrum is the generated by evaluating the best fit spline from the last iteration for the full range of dispersion pixels.

Here is an example of continuum fitting for a section of one order of a stellar spectrum:

.. image:: images/example_continuum_fit.pdf
	:width: 80%
	:alt: Example continuum fit.

The first iteration of the spline fit is heavily biased by the presence of absorption lines, but the final iteration is adequately unbiased by the absorption. It might slightly still underestimate the continuum, but does a relatively good job.

Structure of the continuum extension
------------------------------------

The continuum fit values are appended in an extension to the base spectrum file for each on-sky science observation, which is output by the extraction step (see :ref:`here <target_to_spectrum_file>` for the base spectrum file structure).

The continuum extension is HDU index 4 and named "continuum". Its data entry is the continuum spectra for each order with shape (number of orders, number of dispersion pixels).

The only addition to the primary HDU header is a ``HISTORY`` entry with the date on which continuum fitting is performed.



