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

Below are tables defining the options in each of the main *config* file sections, which govern how the reduction modules are executed.

Each table has three columns: the *config* keyword, a description, and the default value for our base use case with the Tull coudé spectrograph for the echelle setup used by Professor Adam Kraus's group at UT Austin.

**Every module's config block begins with a** ``do_step`` **flag that determines if a module is executed.**

CCD calibration files module
++++++++++++++++++++++++++++

================================= =========================================================================== ==============================
**Config Key**  				  **Description**															  **Default value**
--------------------------------- --------------------------------------------------------------------------- ------------------------------
``bias_image_type``    			  Header image type keyword value for the bias frames						  'zero'
``flat_image_type``    			  Header image type keyword value for the flat frames						  'flat'
``flat_object_names_to_discard``  Header object keyword values to exclude from the creation of the flat field [ 'FF integration time test' ]
``bias_bpm_percentile``    		  Bias flux percentile level above which to mark a bad pixel                  99.9
``flat_field_bpm_limit``          Flat field response value below which to mark a bad pixel                   1.0e-4
================================= =========================================================================== ==============================

Image processing module
+++++++++++++++++++++++

================================= ==================================================================== ==============================
--------------------------------- -------------------------------------------------------------------- ------------------------------
``arc_lamp_image_type``    	      Header image type keyword value for arc lamp spectra                 'comp'
``valid_arc_lamp_object_names``   Header object keyword values to include for arc lamp frames (a list) [ 'thar', 'a' ]
``invalid_science_object_names``  Header object keyword values to exclude as science frames (a list)   [ 'test' ]
``cosmic_subtract``               Flag for whether or not to perform cosmic ray subtraction            True
``cosmic_subtract_niter``         Number of iterations for cosmic ray subtraction                      4
================================= ==================================================================== ==============================

Order tracing module
++++++++++++++++++++

============================== ======================================================================================================================= ==============================
------------------------------ ----------------------------------------------------------------------------------------------------------------------- ------------------------------
``order_xdisp_trace_width``    Cross dispersion pixel height of the slit                                                                               20
``order_start_index``          The index at which to find the starting location of the order traces                                                    -33
``order_center_method``        Algorithm to use for finding the initial locations of the order traces. Must be 'peak_find' or 'gradient_threshold'     'peak_find'
``trace_poly_degree``          Polynomial degree to fit to the trace values                                                                            2
``trace_poly_fit_start_index`` The starting dispersion pixel to use for fitting the trace polynomial                                                   512
``number_of_orders``           Number of orders to trace. If fewer are found, the trace is extended                                                    58
============================== ======================================================================================================================= ==============================

Spectral extraction module
++++++++++++++++++++++++++

=============================== ======================================================================================================================= ==============================
------------------------------- ----------------------------------------------------------------------------------------------------------------------- ------------------------------
``reverse_traced_orders``       Flag for whether or not the order direction needs to be reversed (to match increasing wavelength order)                 True
``order_xdisp_width_extract``   Cross dispersion pixel width of an order for extraction                                                                 -33
``science_extract_type``        Extraction algorithm to use for on-sky science observations. Must be 'optimal_extractin' or 'sum_extraction'            'optimal_extraction'
``lamp_extract_type``           Extraction algorithm to use for arc lamp observations. Must be 'optimal_extractin' or 'sum_extraction'                  'sum_extraction'
``science_background_subtract`` Background subtraction method to use for on-sky science observations. Must be 'fixed' or 'fit'                          'fixed'
``lamp_background_subtract``    Background subtraction method to use for arc lamp observations. Must be 'subtract' or 'none'                            'none'
=============================== ======================================================================================================================= ==============================

Wavelength calibration module
+++++++++++++++++++++++++++++

================================== ======================================================================================================================= ==============================
---------------------------------- ----------------------------------------------------------------------------------------------------------------------- ------------------------------
``min_arc_exp_time``               Minimum exposure time (in seconds) for an arc lamp observation to be included in wavelength calibration                 30
``use_prelim_sol_order_offset``    Flag for whether or not an order-index offset between observation and initial wavelength solution should be found       True
``peak_threshold_mad_method``      Method to use for estimating the arc lamp flux noise for peak finding. Must be 'full_spectrum' or 'chunk_spectrum'      'chunk_spectrum'
``lamp_line_peak_threshold_sigma`` Number of standard deviations above the noise to use for the arc lamp peak finding algorithm                            5
``lamp_line_min_separation_pix``   Minimum separation (in pixels) of two consecutive arc lamp peaks                                                        5
``lamp_line_pix_width_limits``     Range of peak widths (in pixels) to allow in the peak finding algorithm                                                 [ 2, 4 ]
``wave_sol_guess``                 File name of the initial wavelength solution array                                                                      'prelim_wave_sol.npy'
``line_list``                      File name of the reference arc lamp line list wavelengths                                                               'thorium_line_list.csv'
``arc_ref_file``                   File name of the reference arc lamp (photron ThAr) for plotting/order-offset purposes                                   'photron_thar_atlas.csv'
``max_wave_diff_with_list``        Maximum allowed wavelength difference between an observed peak and the reference line list for fit inclusionn           1.0
``wave_cal_poly_order``            Degree of polynomial for the wavelength solution fit                                                                    4
``vel_resid_sigma_reject``         Number of standard deviations to reject velocity residuals in the iterative wavelength solution fitting                 3
================================== ======================================================================================================================= ==============================

Continuum fitting module
++++++++++++++++++++++++

=============================== ======================================================================================================================= ==============================
------------------------------- ----------------------------------------------------------------------------------------------------------------------- ------------------------------
``lower_sigma_reject``          Number of standard deviations to reject flux values below the continuum fit                                             2
``upper_sigma_reject``          Number of standard deviations to reject flux values above the continuum fit                                             5
``num_spectrum_chunks``         Number of chunks to break the spectrum into -- sets the number of spline knots                                          6
``min_knot_spacing``            Minimum knot spacing in wavelength (angstrom)                                                                           15
=============================== ======================================================================================================================= ==============================

Radial velocity module
++++++++++++++++++++++

=============================== ======================================================================================================================= ==============================
------------------------------- ----------------------------------------------------------------------------------------------------------------------- ------------------------------
``orders_to_use_file_name``     File name of the pre-defined orders to use for the broadening function computation                                      'good_orders_tull_coude.csv'
``template_file_name``          File name of the spectrum template to use for broadening function computation                                           'phoenix_t5500_g4.5_m0.0.csv'
``bf_velocity_span``            Velocity span to compute the broadening function over (in km/s)                                                         300
``bf_smooth_res``               Resolution to smooth the broadening function to for fitting (normally instrument resolution)                            60000
``n_bootstrap_samples``         Number of bootstrap samples for brodening function comination to measure the radial velocity                            2500
=============================== ======================================================================================================================= ==============================






