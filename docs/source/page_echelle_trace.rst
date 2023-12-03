Tracing Echelle Orders
======================

.. raw:: html

    <style> .red {color:red} </style>
    <style> .purple {color:#731683;font-weight:bold} </style>

.. role:: red
.. role:: purple

Before extracting 1D spectra from the processed images, we need to define the echelle order traces. The trace defines the center of the spectral order in the cross dispersion direction across the dispersion axis, along which the 1D spectrum can be extracted. In this pipeline, we trace orders using the flat field. The order tracing functionality is hosted in the :py:meth:`trace_echelle <trace_echelle>` module.

The options for this step in the reduction are defined in the ``trace`` section of the main *config* YAML file, which is described in full :ref:`here <target_to_config_description>`.

General outline of order tracing steps
--------------------------------------

An overview of the steps that we implement to trace the echelle orders are:

1. Identify the spectral orders and find their initial center locations at the edge of the image.
2. Trace each spectral order along the dispersion direction from that initial location.
3. Fit the trace cross-dispersion centers as a function of dispersion pixel for each order with a polynomial.
4. Identify bad order traces by fitting the trace polynomial coefficients as a function of order and rejecting outliers.
5. Output the trace file.

We use the flat field to find the order traces because it has high signal, is reliably present every night of observation, and fully fills the slit which provides clear cross-dispersion order edges. In the previous `Python 2 version of this reduction pipeline <https://github.com/dkrolikowski/coudereduction>`_, we found the initial trace locations with the flat and then the traced the full orders using a combination of bright science frames taken throughout the night. This was unreliable due to the variably quality of science frames taken on a night (and would not be possible on nights without science frames, but which might still have arc lamp frames that are of interest for wavelength calibration analysis).

Getting the order traces
------------------------

The first thing we need to do is find the pixels that trace the cross-dispersion centers of each order of the echellogram. We do this by first identifying the orders at the edge of the image, and then following each order "block image" across the dispersion direction to get the full order traces. 

Finding the initial trace locations
+++++++++++++++++++++++++++++++++++

We identify the initial order locations using a slice along the cross-dispersion axis of the detector image ("y"-axis) at its edge, just before the overscan region. 

Below are plots to highlight these steps:

.. image:: images/trace_start.pdf
	:width: 90%
	:alt: Plots to illustrate the identifcation of order and initial starting trace locations.

The left panel shows a portion of the flat field image. The blue vertical line at the right edge of the image shows the dispersion pixel at which the vertical slice of the flat field is taken to find the initial location of the orders. This is right inside the overscan region. We start at the right edge of the detector because the orders curve downwards as they move to the left in the image, meaning that only orders present at the right edge are fully on the detector. The image highlights how the flat lamp fills the entire slit, producing clear edges at the cross-dispersion top and bottom of each order.

The right panel shows the flat field values along the cross-dispersion slice at the right edge of the detector. The orders are shown as the rectangular peaks across the slice -- these are what we want to find the centers of as the starting points for the order traces.

To find the order peaks, we use the gradient of the slice of the flat field. At the sharp order edge, the gradient of the flat field peaks sharply as the order edge rises (and similarly at the other edge of the order the gradient negatively peaks as the order decreases). Thus, we can look for the positive peaks of the gradient to denote the start of an order, and then find the center of each identified order.

There are two different algorithms implemented to find the flat field slice gradient peaks in the :py:func:`find_order_centers_along_slice <trace_echelle.find_order_centers_along_slice>` function:

1. Direct thresholding based on an estimate of the noise of the flat slice's gradient to identify peaks above a certain level. The noise is estimated using the median absolute deviation of the gradient of the flat slice at the bottom of the detector, where the response is low and no orders are present above the noise.
2. Using the ``scipy.signal`` ``find_peaks`` algorithm (`documentation here <https://docs.scipy.org/doc/scipy/reference/generated/scipy.signal.find_peaks.html>`_) with constraints placed on the width of the peaks. The range of allowable peak widths is hard coded into the function and was determined by directly inspecting a flat field from the Tull coudé spectrograph. We also set any values of the flat slice gradient below its median to the median, as the ``find_peaks`` algorithm has issues when the negative gradient peaks are present (which is okay because we only need the starting edge of the orders)

The algorithm that is used is defined in the *config* file.

The pipeline default is to use the ``scipy`` function. This is more reliable than the direct thresholding, which is hard to tune for the varying flat field response resulting in significantly different gradient peak values for each order. The ``scipy`` function does require hard coding of the peak width constraint, although in the future that can be changed to a *config* file option. The ``scipy`` function much more reliably finds the flat slice gradient peaks, and also does a better job of finding more orders towards the bottom of the detector where the signal signficantly degrades. 

**Important note**: The flat slice gradient peak finding identifies what is roughly the starting edge of an order. However, we want to identify the centers of the orders. In the *config* file we define the cross-dispersion pixel height of the slit, and add half of that value to the peak finding output to translate them to order centers. We then re-center the order locations by identifying the edges of the order flat slice as being where the values are 70% of the maximum, and then adopt the halfway point as the order center. This re-centering is done with the :py:func:`recenter_order_trace <trace_echelle.recenter_order_trace>` function.

:purple:`Diagnostic plot:` a plot like the above figure's right panel showing the flat field slice and the centers of the identified orders is output in the ``trace`` subdirectory.

Tracing the order across the detector
+++++++++++++++++++++++++++++++++++++

With the orders identified and their centers at the edge of the detector measured, we can trace the centers of each order across the full detector to get the full traces. This is done with the :py:func:`find_full_trace <trace_echelle.find_full_trace>` function.

We start with the initial centers found, then move a pixel to the left and re-center the trace as is done at the end of the initial trace location finding step. We iterate this for every dispersion pixel and for each order that is identified. This results in an array of trace centers with shape (number of orders identified, number of dispersion pixels). For the Tull coudé spectrograph the latter is 2048 pixels.

It is possible for some orders to be poorly traced. This is particularly true towards the bottom of the detector where signal significantly degrades. There is also an artifact on the detector, called the "picket fence", that imprints an emission like fringe over some orders. This "picket fence" fringing can cause those orders' traces to be poor. This is fixed in the fitting of the trace.

:purple:`Diagnostic plot:` a plot showing the flat field with the full found trace overplotted as points is output in the ``trace`` subdirectory. It is a multi-page figure, with one page showing the full detector, and two pages to zoom in on each the top and bottom half to better see the trace. Generated with the :py:func:`plot_trace_on_image <trace_echelle.plot_trace_on_image>` function.

Fitting the trace
-----------------

We then fit the full found trace across the detector for each order with a polynomial to enfore smoothness. For each order we fit the cross-dispersion pixel center vs. dispersion pixel. The degree of the polynomial is defined in the *config* file, and our default is 2nd-order. There is also an option in the *config* file to set the starting dispersion pixel to fit the polynomial to. This was introduced because sometimes the trace at the left edge of the detector wanders off due to low signal and biases the fit, despite the rest of the trace being identified well. For our default use, this is set to exclude the first fourth of the order.

:purple:`Diagnostic plot:` a similar plot to that for the full trace finding is output to the ``trace`` subdirectory, with the polynomial fit to each order's trace overplotted. Generated with the :py:func:`plot_trace_on_image <trace_echelle.plot_trace_on_image>` function.

Vetting quality of order trace fits
+++++++++++++++++++++++++++++++++++

As stated above, some orders can have poorly defined traces. This would result in bad fits to their traces, and would produce spurious 1D extracted spectra for those orders. We vet the quality of the order traces by looking at the trace fit polynomial coefficients across orders, which should be smooth. We fit each of the polynomial coefficients with their own 2nd-order polynomial as a function of order, and any orders with coeffiecients greater than :math:`10\sigma` (as calculated by the median absolute deviation of the residuals) away from the 2nd-order polynomial fit are marked as bad orders. An order only needs one bad coefficient to be marked as a bad order.

These bad orders have their trace fit polynomial coefficients replaced with the value of the "hyper-fit" at that order value (zero-indexed). If there are bad order traces, the orders are marked in a ``HISTORY`` entry of the trace FITS file. 

Here is an example showing a fit to the 0th order coefficient (intercept) of the trace fit polynomial as a function of order. The bad orders, of which there are 4, are shown as the brown points. The 3 bad orders at end are due to low signal of the flat field and the bad order in the middle is due to the "picket fence" contamination.

.. image:: images/example_trace_coeff_hyperfit.pdf
	:width: 90%
	:alt: Example hyper-fit to the 0th order coefficient of the trace fits.

The initial fitting of the trace and subsequent vetting of the polynomial fits is done with the :py:func:`fit_full_trace <trace_echelle.fit_full_trace>` function.

:purple:`Diagnostic plot:` a multi-page plot showing a similar figure as above for each coefficient on its own page is output in the ``trace`` subdirectory.

Final adopted trace fit
+++++++++++++++++++++++

The final adopted fit trace values are then the values of the full trace center vs. dispersion pixel polynomial fits. These polynomial coefficients are only replaced with the result of the "hyper-fit" as a function of order for the identified bad orders, otherwise the order's individual fit is used.

The *config* file also sets the number of orders to trace, so the user has control over which orders are extraced. As our default, this value is 58. If fewer orders are found in the beginning of the trace step than requested, the trace polynomial coefficient fits as a function of order are used to extend the fit trace (in this case, these extended orders are added to the list of bad orders in the ``HISTORY`` entry).

Here is an example plot showing a zoom-in of the final adopted fit trace plotted on top of the full found trace and flat field (for the bottom half of the detector):

.. image:: images/example_final_trace.pdf
	:width: 70%
	:alt: Example final adopted trace.

The pink points are the data-found trace points and the lines are the polynomial fits to them. The red solid lines denote orders that have good trace fits and the dashed teal lines denote bad orders whose trace fits have been replaced with coefficients from the "hyper-fit" to the coefficients vs. order. The bad order near the top of the image, a result of the "picket fence", highlights how the vetting of the trace fits is necessary to provide a good order trace.

:purple:`Diagnostic plot:` a plot similar to the above figure -- showing the trace data, fits, and flat field -- is output to the ``trace`` subdirectory, but with multiple pages like the diagnostic plot for the initial found and fit traces to show the full detector and zoom-ins of the top and bottom half. Generated with the :py:func:`plot_trace_on_image <trace_echelle.plot_trace_on_image>` function.

Trace FITS file structure
-------------------------

The order trace is saved to a file called "trace.fits" in the ``trace`` subdirectory.

It is a FITS file with 3 extensions:

- Index 0: The primary extension that only contains header information.
- Index 1: Extension named "data trace" with the order centers traced along the dispersion direction. It has shape (number of orders, number of dispersion pixels).
- Index 1: Extension named "fitted trace" with the order centers computed by the polynomial fit to "data trace". It has shape (number of orders, number of dispersion pixels).

The primary extension has additional keywords relevent to the trace finding:

============ ===================================================================================================================================================================================
**Keyword**  **Description**
------------ -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
``NORDERS``  the number of orders traced.
``POLYDEG``  the polynomial degree used for fitting the trace.
``HISTORY``  multiple entries with details of various other options and parameters used (e.g. intial order center finding algorithm used, bad orders whose trace fit parameters were corrected.)
============ ===================================================================================================================================================================================


