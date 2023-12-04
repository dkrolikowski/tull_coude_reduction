Wavelength Calibration
======================

.. role:: purple

The last step of the "reduction" part of the pipeline to produce science-ready spectra is the wavelength calibration. We do this using observations of arc lamps, which are ThAr lamps in our default case with the Tull coudé spectrograph. We use a reference line list to map the pixel centroids of identified emission lines in the observed arc lamp spectra to their corresponding wavelengths in the line list. Then a polynomial is fit to produce a pixel -> wavelength solution for the entire order, and applied to every science observation.

The :py:meth:`wavelength_solve_and_calibrate <modules.wavelength_solve_and_calibrate>` module contains the functions for wavelength calibration. The options for this step are defined in the ``wavecal`` section of the main *config* YAML file, which is described in full :ref:`here <target_to_config_description>`.

.. warning::

	As the pipeline currently stands, the wavelength calibration module is by far the least flexible.

	Creating a bespoke wavelength solution on-the-fly for any arbitrary echellograph spectrograph set up is a very difficult problem. We circumvent this by necessitating an initial wavelength solution for every order as input. Thus, at least once a wavelength solution needs to be generated manually (or at least separate from the pipeline). As a result, this wavelength calibration algorithm is optimized for a lot of data taken with the same echelle set up.

	Working towards a more flexible module is an area of active development.

	There is also a script to interactively define a wavelength solution in the original Python 2 version of the pipeline that will be updated for use in making new initial wavelength solution guesses for different echelle set ups.

Wavelength solution from arc lamps
----------------------------------

The first part of wavelength calibration is to measure a wavelength solution mapping pixel to wavelength using a reference source -- in this case an arc lamp spectrum. This step is run with the :py:func:`wavelength_solution <modules.wavelength_solve_and_calibrate.wavelength_solution>` function.

Reference files
+++++++++++++++

There is a set of reference files that is necessary for the wavelength calibration step. These are defined in the *config* file and are kept in the ``data`` directory in the ``tull_coude_reduction`` code repository.

- A line list with catalogue wavelengths of emission lines for the arc lamp that is observed. In the default Tull coudé case this is the Redman thorium line list for the ThAr lamp. This is saved as a single column CSV file with a column named "wavelength".
- An initial wavelength solution guess. This provides an initial wavelength solution that is refined for any drift, such as those induced by echelle setup shifts. This is saved as a .npy file that is a 2D ndarray object, with shape (number of orders, number of dispersion pixels). It is in essence a wavelength calibration as output by this module of the pipeline.
- A reference wavelength-calibrated ThAr spectrum. This is an observation of a Photron ThAr lamp (`from KPNO atlas <https://noirlab.edu/science/data-services/other/spectral-atlas>`_), and is a CSV with columns "wavelength" and "flux". This is used to measured an order offset between the initial wavelength solution guess and the night's observations, and in diagnostic plotting.

Preparation of the arc lamp observations
++++++++++++++++++++++++++++++++++++++++

We first need to identify the arc lamp observations to use for measuring the wavelength solution. The ``image_process`` *config* section sets which observations are arc lamps, using the ``imagetyp`` FITS header keyword (default here is "comp") and a list of ``object`` FITS header values that denote valid arc lamp spectra (in our case, this excludes arc lamp frames taken for spectrograph focusing). We also include an option in the ``wavecal`` config section that sets a minimum exposure time for arc lamp observations to include in calibration. The default value is 30 seconds, which excludes shorter exposure times that have much too little signal for accurate calibration, particularly in the blue.

With the arc lamp observations chosen, we can go frame by frame and derive a wavelength solution for each one.

Initial wavelength solution order offset
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The pipeline has the option, set in the *config*, to determine if there is an order-index offset between the initial wavelength solution guess and the night's observations. One reason this would be needed is if the initial solution has more orders than is extracted in the pipeline -- this is true for the default setup, which had a initial solution created with order padding. The offset resets the 0th order index of the initial solution. Thus, it is not allowed to be negative or call for an order beyond the size of the initial solution. This offset is determined with the :py:func:`order_offset_with_wave_sol_guess <modules.wavelength_solve_and_calibrate.order_offset_with_wave_sol_guess>` function.

To measure the offset, we take an order of the observed arc lamp spectrum and assign it wavelength solutions from the initial guess with a range of offsets from its order index. We calculate the residual between the observed spectrum and the reference Photron ThAr spectrum (both median-normalized) and adopt the order offset with the minimum residual. This is repeated for multiple orders and the most common order offset is used moving forward (if there are multiple calculated). By using multiple orders, we reduce the effect of low signal.

In our default use case, the initial solution has 59 orders but we extract 58. The order offset is 1, meaning that order index 0 in the observed spectrum corresponds to order index 1 in the initial solution.

Identifying emission lines
++++++++++++++++++++++++++

With the order offset for the initial solution determined, we can now wavelength calibrate the arc lamp spectrum. To do this, we need to go order by order and identify emission lines that will be matched to the reference line list to map pixel to wavelength.

The function :py:func:`find_arc_lamp_line_pixel_centers <modules.wavelength_solve_and_calibrate.find_arc_lamp_line_pixel_centers>` takes in an order of the extracted spectrum and finds emission lines as peaks. We normalize the extracted spectrum by its median to make it easier to determine which peaks are significant. 
Similar to the trace finding, we use the ``scipy.signal`` ``find_peaks`` algorithm to identify emission lines (`documentation here <https://docs.scipy.org/doc/scipy/reference/generated/scipy.signal.find_peaks.html>`_).

Three different constraints in the ``find_peaks`` algorithm are used to find the peaks: distance, width, and prominence:

- **Distance**: This defines the minimum separation between consecutive peaks in pixels and is set in the *config* file. This ensures that blended lines are not included in the wavelength calibration.
- **Width**: This sets the limits on the pixel width of the peaks and is set in the *config* file. An upper constraint on the width also helps to exclude blended lines, and the lower constraint helps to reduce the number of oxide lines in a contaminated ThAr lamp that are included in the peak list.
- **Prominence**: This sets the minimum prominence (a signal-to-noise) for a peak to be included. By normalizing the spectrum by the median, this is essentially the height of the line. We estimate the noise of the spectrum using the median absolute deviation of the flux, and set a number of standard deviations above that noise in the *config* as the minimum prominence. We use the :py:func:`get_flux_mad_from_spectral_chunks <modules.wavelength_solve_and_calibrate.get_flux_mad_from_spectral_chunks>` function to estimate the order's noise using chunks of the spectrum to reduce the influence of concentrated emission lines (e.g. oxide bands).

Here is an example of a portion of a ThAr spectrum with the identified peaks marked with vertical lines:

.. image:: images/example_thar_spectrum.pdf
	:width: 40%
	:alt: Example ThAr spectrum with found peaks.

We then go through each of the identified peaks and fit it with a Gaussian to get the fractional pixel centroid. These pixel centroids are the candidate arc lamp emission lines that will be matched to the line list and used for wavelength calibration.

Fitting a wavelength solution
+++++++++++++++++++++++++++++

With a list of pixel centroids for candidate arc lamp emission lines, we use the :py:func:`fit_wavelength_solution <modules.wavelength_solve_and_calibrate.fit_wavelength_solution>` function to match the emission lines to the reference line list and fit a polynomial wavelength solution.

We use the initial wavelength solution to map the pixel centroids to wavelength centroids and calculate the difference between each peak's wavelength centroid and the reference line list. Peaks that have wavelength centroids close enough to a line in the reference line list are marked as true arc lamp emission lines and assigned the catalogue wavelength value. The *config* file sets the maximum wavelength difference allowed between the centroid and reference line list to still be called a match.

A polynomial is fit to the set of peak pixel centroids and their corresponding line list wavelengths. The degree of the polynomial is set in the *config*. Then, velocity residuals are calculated between the wavelength polynomial solution and line list wavelength. Lines with a velocity residuals larger than a *config* defined number of standard deviations (as measured with the median absolute deviation of the residuals) are rejected. This fitting is iterated until no lines are rejected or the number of remaining lines would be smaller than the degree of the polynomial + 1.

Here are example plots showing one spectral order's velocity residuals for two iterations of fitting (the first and last iteration).

.. image:: images/example_wavefit_residuals.pdf
	:width: 95%
	:alt: Velocity residuals for an example ThAr spectral order.

The sinusoidal residual in the first iteration (left panel) is the result of an offset between the initial wavelength solution and the actual wavelength solution for the night's echelle setup. The last iteration (right panel) has relatively low amplitude residuals without structure, highlighting the quality of the solution after the iterative fitting routine.

The wavelength solution polynomial is then evaluated for each order, and is output in the arc lamp observation's spectrum file as its wavelength solution. See the wavelength extension structure section below for more details about the file output.

:purple:`Diagnostic plot:` One diagnostic plot for the wavelength fitting is a multi-page PDF showing the velocity residuals for each iteration, like the plots shown above. Each page is one iteration, with the order reversed: the first page shows the final adopted iteration and so on. The plot is saved in the ``wavecal/fit_residuals`` subdirectory. Each order has its own plot. Created with the :py:func:`plot_wavelength_fit_iteration_residuals <modules.wavelength_solve_and_calibrate.plot_wavelength_fit_iteration_residuals>` function.

:purple:`Diagnostic plot:` Another diagnostic plot: a multi-page PDF showing the observed arc lamp spectrum with the reference spectrum overplotted. The lines included in the wavelength solution fit are marked with vertical dashed lines. Each page is one iteration, with the order reversed: the first page shows the final adopted iteration and so on. The plot is saved in the ``wavecal/fit_residuals`` subdirectory. Each order has its own plot. Created with the :py:func:`plot_wavelength_fit_iteration_spectra <modules.wavelength_solve_and_calibrate.plot_wavelength_fit_iteration_spectra>` function.

:purple:`Diagnostic plot:` Another diagnostic plot: a similar set of figures to the previous plot except with zoom-ins of spectral windows using the adopted wavelength solution to provide a fine-grained look at the quality of the fit. Each window covers 10 angstrom and a maximum of 6 windows are included on a page. Each plot is multi-page for one order. Created with the :py:func:`plot_spectra_zoom_windows <modules.wavelength_solve_and_calibrate.plot_spectra_zoom_windows>` function.

Wavelength calibration of science frames
----------------------------------------

The wavelength solution fitting provides wavelength calibration for the arc lamp observations. However, we need to wavelength calibrate the rest of the science observations taken throughout the night. 

The wavelength calibration step provides pairs of observation times and wavelength solutions for the arc lamp spectra, which are typically taken at the start and end of the night (and sometimes in the middle of the night). We linearly interpolate the arc lamp solutions to each science frame's observation time to provide a wavelength solution. If a science observation is not bracketed by arc lamp observations, we then just adopt the arc lamp solution that is closest in time to the science observation.

The interpolation is done with the :py:func:`interpolate_wavelength_solution <modules.wavelength_solve_and_calibrate.interpolate_wavelength_solution>` function.

Structure of the wavelength extension
-------------------------------------

Each science observation has a wavelength extension appended to its spectrum file, which is output by the extraction step (see :ref:`here <target_to_spectrum_file>` for the base spectrum file structure).

The wavelength extension is HDU index 3 and named "wavelength". Its data entry is the wavelength solution with shape number of orders, number of dispersion pixels.

The primary HDU of the spectrum file has new keywords added to it:

============ =================================================================================================================================
**Keyword**  **Description**
------------ ---------------------------------------------------------------------------------------------------------------------------------
``WAVPOLYD`` the polynomial degree of the wavelength solution polynomial.
``WAVTYPE``  the interpolation used for generating the wavelength solution (either "linterp" for linear interpolation or "closest" if the closest calibration is adopted). **Only for science frames not used for arc lamp calibration (on-sky observations and arc lamps with too little exposure time).**
``HISTORY``  one entry containing the date that the wavelength calibration was done.
============ =================================================================================================================================

Note that the headers are slightly different between the arc lamp observations used for calibration and the rest of the science frames. Only the latter have the ``wavtype`` keyword, because the arc lamp observation used for calibrations do not require interpolation to generate a solution (because they themselves are fit!)






