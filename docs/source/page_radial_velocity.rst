Measuring stellar radial velocities
===================================

.. role:: purple
.. role:: blue

At this point in the pipeline, we can use the continuum measurement to normalize the 1D extracted spectra to measure radial velocities (RV). We use the ``saphires`` package (`documentation here <https://saphires.readthedocs.io/en/latest/intro.html>`_) to compute the broadening function of the observed stellar spectrum using a model template spectrum. With the broadening function we can measure an RV, with precisions below a km/s for many targets (and even to below 100 m/s for the brightest targets).

See the ``saphires`` documentation, and the Rucinski references mentioned there, for more information about broadening functions -- we won't discuss the details of broadening functions themselves below.

The :py:meth:`radial_velocity <modules.radial_velocity>` module contains the functions for measuring radial velocities from the extracted, continuum normalized 1D spectra. The options for this step are defined in the ``radial_velocity`` section of the main *config* YAML file, which is described in full :ref:`here <target_to_config_description>`.

Preparation
-----------

Barycentric velocity correction
+++++++++++++++++++++++++++++++

The velocity shift of the observed spectrum is the relative velocity between the observer and the observed star. To get the actual RV of the star, the motion of the earth needs to be removed from the measured velocity shift. This is called a barycentric velocity correction. In this pipeline we calculate the barycentric velocity correction using the ``barycorrpy`` package, which is described in its `documentation here <https://github.com/shbhuk/barycorrpy/wiki>`_.

The barycentric velocity correction calculation requires an observation time, the on-sky coordinates of the target, and observatory location. We take the observation time from the primary extension header. If there is exposure meter data input into the header, we use the ``EMFWMTIM`` keyword value, which is the exposure meter flux weighted midtime of the observation. If this keyword is not in the header, then the plain midtime of the observation is used. The on-sky coordinates are taken from the ``RA`` and ``DEC`` primary extension header values. The observatory location is set to use the preset "McDonald" observatory values from ``barycorrpy``.

The barycentric velocity correction in km/s is added to the spectrum FITS file primary extension header under the keyword ``BARYCORR``. 

.. warning::

	It is possible that the header right ascension and declination values aren't correct. We haven't fully evaluated how common a mistake is, or by how much the values might be wrong. This would affect the accuracy of the barycentric velocity correction, although likely not to any effect greater than the typical Tull spectrograph RV precision. We will work on estimating the possible range of systematic barycentric velocity correction error that could be introduced from this.

Orders to use
+++++++++++++

The spectral orders for which we compute broadening functions are currently defined in a data file in the package repository. The file name is defined in the *config* file as ``orders_to_use_file_name``. This file is a CSV with four columns: the order index ("order"), starting wavelength of order ("wave_start"), ending wavelength of order ("wave_end"), and a flag ("flag"). The flag can have values of either "start" or "end" and denotes orders for which the full extent shouldn't be used. For example, if the flag is set to "start", the order will only be used from its "wave_start" value onwards. These flags are used in orders with either significant telluric bands or very strong stellar lines (e.g. H-alpha) to avoid.

In the future the choosing of which orders to use will be made automatic. There still should be some constraints: a starting and ending wavelength (to avoid low signal orders for example), a telluric model to avoid strong telluric bands, and strong stellar features that should be excluded. However, this can be implemented in a manner that is agnostic to the echelle setup so that it is flexible for any use of the Tull spectrograph.

Broadening function: computation and analysis
---------------------------------------------

Computation steps: process the template, detail what is currently used. The compute broadening functions for each order.

Measuring RVs: bootstrap sample, combining different orders, smoothing the broadening function. Fit with a gaussian, adopt the median and scaled-MAD as the value and error. Errors may be under-estimated.

:purple:`Diagnostic plot:` Two panel. Left panel is the combined broadening function, a gaussian fit, and two vertical lines showing the bootstrap RV + 3sigma range and the barycentric velocity. Right panel is the bootstrap RV distribution with vertical line denoting the median and shaded region the 1 sigma range.

Known issues
++++++++++++

Could be sky contamination that affects the BF, see a peak at the barycentric velocity. So be careful!

.. warning::

	The sky contamination that seems to be a problem.

Also errors haven't been verified, may be underestimated. Also incorporate rotational broadening measurement.

Structure of radial velocity extension
--------------------------------------

The radial velocity results are appended to the base spectrum file for each on-science observation. This file is initially output by the extraction step (see :ref:`here <target_to_spectrum_file>` for the base spectrum file structure).

The radial velocity extension is HDU index 5 and named "radial velocity". Its data entry is the order-by-order smoothed broadening function with shape (number of orders, number of broadening function velocity points). Orders for which a broadening function is not computed are given an array of ``nan``. The velocity array is the same for all order broadening functions and can be recreated with information in the extension header, which provides the first value in the velocity array, the velocity spacing of the array, and the number of points. :blue:`Note that the broadening function has not been barycentric velocity corrected.` 

The RV value and error, calculated with the bootstrap samples, are output in both the primary extension and radial velocity extension headers. :blue:`The RV value HAS been corrected for the barycentric velocity.`

The "radial velocity" extension has these keywords related to the BF computation and radial velocity measurement:

============ =================================================================================================================================
**Keyword**  **Description**
------------ ---------------------------------------------------------------------------------------------------------------------------------
``RVBF``     the measured radial velocity in km/s (median value of bootstrap sampling)
``ERVBF``    the measured radial velocity error in km/s (standard deviation of bootstrap sampling)
``NRVBOOT``  the number of bootstrap samples used for RV measurement
``VELSTART`` the starting value for the broadening function velocity array in km/s
``VELSTEP``  the spacing of the brodening function velocity array in km/s
``NVELPTS``  the size of the broadening function velocity array
============ =================================================================================================================================

Note that ``RVBF`` and ``ERVBF`` are also added to the primary extension's header. The only addition to the primary HDU header ``HISTORY`` is an entry with the date on which the RV processing occurred.


