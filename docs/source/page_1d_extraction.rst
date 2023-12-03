Extraction of 1D Spectra
========================

.. role:: purple

Now that location of the echelle orders have been found and traced, we can extract 1D spectra from the 2D spectrum image of each order. This is done by summing up the flux in the cross-dispersion slice of the order flux image at every dispersion pixel, collapsing the 2D spectrum into one dimension. The result of this step is a series of 1D spectra, flux and flux error vs. dispersion pixel, for each order that was traced in the previous step.

The :py:meth:`extract_spectrum <extract_spectrum>` module contains the functions for extracting the spectra. The options for this step are defined in the ``extraction`` section of the main *config* YAML file, which is described in full :ref:`here <target_to_config_description>`.

Extraction algorithms
---------------------

Two extraction algorithms are current included for use in the pipeline: sum extraction and optimal extraction. The pipeline is set up so that different extraction algorithms can be used for different types of science frames. In our default case, this differentiates arc lamp observations (emission line spectra) and on-sky science observations (stellar absorption spectra).

The first step, regardless of extraction method, is to cut the flux and error 2D image blocks for each order from the full processed detector image. We do this with the :py:func:`get_order_image_block <extract_spectrum.get_order_image_block>` function. For each order the function goes pixel by pixel in the dispersion direction, rounds the fit trace value at that dispersion pixel, and defines the cross-dispersion slice of the 2D spectrum with a *config* defined cross-dispersion order pixel width (default: 16 pixels).

With the order flux and error image blocks, we then can go pixel by pixel and collapse the 2D spectrum into 1D spectra. Below we describe the available extraction algorithms for this step.

Optimal extraction
++++++++++++++++++

The function :py:func:`optimal_extraction <extract_spectrum.optimal_extraction>` is an implementation of the Horne 1986 optimal extraction algorithm.
We will not discuss the details behind the algorithm, only the broad overview of steps taken and the specifics/differences of our implementation. See the original `Horne 1986 <https://ui.adsabs.harvard.edu/abs/1986PASP...98..609H/abstract>`_ paper for the details of optimal extraction.

The general idea is to assume a spatial profile of the flux in the cross-dispersion direction, fit each cross-dispersion slice of the 2D order flux image block, and use the normalized spatial profile to calculate the summed flux in each dispersion pixel. Here we assume a Gaussian spatial profile, which is defined in the :py:func:`gaussian_1d <reduction_utils.gaussian_1d>` function. 

We loop through each dispersion pixel and fit the cross-dispersion slice of the order flux image block with a Gaussian. There are two currently available options for setting the background flux value, which is set in the *config* file: ``fit`` or ``fixed``. See the subsection below for more information about the background subtraction options, but each of the options available provides an offset flux level for the flux slice Gaussian. The Gaussian fit is performed with ``scipy.optimize.curve_fit``, and we include the cross-dispersion slice of the order error image block in the fitting routine. We exclude ``nans`` (bad pixels) from the fit, and if there are fewer than 4 non-``nan`` values to fit we skip the slice and the extracted spectrum is set to ``nan``. 

Here is an example of the Gaussian fit to a cross-dispersion slice of an order for an observation of a star:

.. image:: images/example_xdisp_slice_gauss_fit.pdf
	:width: 50%
	:alt: Example Gaussian fit to x-disp pixel slice.

A Gaussian is a reasonable approximation of the spatial profile here.

We then turn the Gaussian fit into a spatial profile by subtracting the background value and normalizing by the sum of the Gaussian fit evaluated at each cross-dispersion pixel (the red star markers in the above plot).

To improve the signal of the extracted spectrum, we can combine the spatial profile parameter information from all dispersion pixels of an order and enforce smoothness across the order. This also allows us to define an accurate spatial profile at dispersion pixels with poor spatial profile fits (or significant outlier flux values along their cross-dispersion slices, for example).

To do this, we fit the spatial profile parameters as a function of dispersion pixel with a polynomial, with sigma rejection iteration. The fitting routine is in the :py:func:`polynomial_fit_sigma_reject <reduction_utils.polynomial_fit_sigma_reject>` function. We use a 3rd degree polynomial with 1 iteration of 5-:math:`\sigma` rejection. We fit the spatial profile Gaussian amplitude and standard deviation parameters, and the background value if the *config* option is ``fit``. The fitting function also takes as input valid ranges of values for each parameter to exclude extreme/unphysical outliers from the fitting (e.g. amplitude values outside of the range 0 to 1).

Here are examples of full-order fits to the spatial profile Gaussian amplitude and standard deviation paramters for one order of an observation of a star:

.. image:: images/example_spatial_prof_par_fits.pdf
	:width: 90%
	:alt: Example fit to spatial profile parameters across dispersion pixel.

Underneath pixel-to-pixel noise, the fit parameter values **do** exhibit a smooth relation across dispersion pixel. 

We again loop through each dispersion pixel to actually extract the 1D spectrum. For each slice, we evaluate the full-order spatial profile paramter polynomial fit at that dispersion pixel value and generate the Gaussian spatial profile. We then use the optimal extraction formula from Horne 1986 to get the 1D extracted flux value, :math:`f`, where :math:`i` is the cross-dispersion pixel along the slice, :math:`b` is the background flux value, :math:`p` is the spatial profile, and :math:`\sigma_f` is the flux error slice:

.. math::
	
	f = \sum_{i} \frac{ ( f_{i} - b ) * p_{i}}{\sigma_{f,i} ^ 2} * \left ( \sum_{i} \frac{p_i^2}{\sigma_{f,i} ^ 2} \right )^{-1}

We also follow the Horne 1986 equations (Equations 9 and 13) for calculating the extracted error spectrum.

.. warning::

	The extracted flux error spectra have not been vetted in great detail. From a cursory investigation, it appears as though the extracted errors may be overestimated.

	There will be future work to quantify this overestimation and track down where the error calculation may have issues.

Background subtraction
^^^^^^^^^^^^^^^^^^^^^^

We have two options for background subtraction with optimal extraction, set in the *config* file: ``fit`` and ``fixed``.

For ``fit``, the background flux value for a cross-dispersion slice is left free to be fit as the Gaussian's offset. 

For ``fixed``, we fit polynomials to dispersion slices at the top and bottom of the order flux image block to estimate the background values. We use :py:func:`polynomial_fit_sigma_reject <reduction_utils.polynomial_fit_sigma_reject>` to fit 2nd order polynomials with 1 round of 5-:math:`\sigma` rejection (to exclude outlier flux values, such as sky emission lines). Then, the background value for a cross-dispersion slice is taken to be the average value of these two background polynomial fits at that dispersion pixel.

The default is to use the ``fixed`` background subtraction. This prevents spurious extracted spectrum pixels where there are sky emission lines, which fill the slit and have their background values fit far too high. In the future, this ``fixed`` background subtraction will essentially occur earlier in the pipeline, with background/scattered light subtraction in the image processing step.

Sum extraction
++++++++++++++

The sum extraction algorithm is relatively simple in its implementation: we loop through each dispersion pixel, subtract a background value from the cross-dispersion flux image slice (if the *config* is set to), and then sum the cross-dispersion pixel fluxes. The extracted flux error is taken from the sum of the squared error image slice, with the background flux value added in if a background is subtracted.

We currently allow for two background subtraction options for sum extraction, which is set in the *config*: ``none`` or ``subtract``. If ``none`` is chosen, no background is subtracted. If ``subtract`` is chosen, then the background flux value is estimated as in the ``fixed`` option for optimal extraction (see above) and subtracted.

Sum extraction without background subtraction is the default for arc lamp spectra.

.. note::

	The benefits of optimal extraction are most clear for faint targets with relatively low signal to noise. While sum extraction is very quick, the computation cost of this optimal extraction algorithm isn't significant enough to warrant turning it off for brighter science targets. Therefore, we recommend keeping optimal extraction as the default for all on-sky observation extraction.

Improvements to make
++++++++++++++++++++

There are improvements that can be made to the above extraction algorithms to produdce even higher fidelity 1D spectra:

- We need to better handle additional bad pixels -- those that get through the cosmic ray subtraction or bad pixel mask. This could be done in a residual rejection scheme after fitting the spatial profile in the optimal extraction.
- Background subtraction is still somewhat unreliable, although the ``fixed`` option for optimal extraction circumvents issues present in the previous version of this pipeline. However, a true scattered light/background subtraction needs to be implemented in the image processing step, which will ease the handling of the background here.

.. _target_to_spectrum_file:

Extracted spectrum file structure
---------------------------------

Now that we have extracted spectra for every traced order we need to output it to a FITS file. The extracted spectra files are saved in the ``spectrum_files`` subdirectory, with file name formatted as "tullcoude_YYYYMMDDTHHMMSS_spectrum.fits", where YYYYMMDDTHHMMSS is the unique observation file token.

These extracted spectra files are FITS files with 3 extensions:

- Index 0: The primary extension only contains header information, which includes a copy of the header from the processed image file.
- Index 1: Extension named "extracted flux" whose data entry is the extracted flux spectra with shape (number of orders, number of dispersion pixels).
- Index 2: Extension named "extracted flux error" whose data entry is the extracted error spectra with shape (number of orders, number of dispersion pixels).

The primary extension has additional keywords relevent to the extraction:

============ =================================================================================================================================
**Keyword**  **Description**
------------ ---------------------------------------------------------------------------------------------------------------------------------
``NORDERS``  the number of orders extracted.
``EXTRACT``  the extraction algorithm used (either "sum_extraction" or "optimal_extraction")
``BGSUB``    the background subtraction option used (either "none" or "subtract" for sum extraction, "fixed" or "fit" for optimal extraction).
``ORDWIDTH`` the cross dispersion order width set in the *config* file.
``HISTORY``  one entry containing the date that the extraction processing was done.
============ =================================================================================================================================

Importance of extracted spectrum files
++++++++++++++++++++++++++++++++++++++

This spectrum file will be the base file that subsequent pipeline steps add to, for now including:

- The wavelength calibration (for all science frames)
- Continuum fitting result (for on-sky observations, if the step is marked in the *config* to do)
- Radial velocity measurement (for on-sky observations, if the step is marked in the *config* to do)







