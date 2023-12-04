Pipeline Overview and Structure
===============================

.. role:: purple
.. role:: blue

Overview of the pipeline steps
------------------------------

This pipeline has two main components: the "reduction" steps to produce wavelength-calibrated 1D extracted spectra and "analysis" steps to measure useful quantities from the spectra for research purposes.

The general steps of the reduction pipeline are:

1. Generate bias, flat field, and bad pixel mask CCD calibration files, and then image process the night's raw CCD images.
2. Trace the location of the echelle orders.
3. Extract 1D spectra along the echelle order traces for science observations.
4. Measure a wavelength solution using reference lamp observations (e.g. ThAr) and calibrate all science observations.

After reduction, the following analysis steps are included:

1. Fit the continuum of stellar spectra.
2. Measure stellar radial velocities using broadening functions.

Pipeline structure: how to run it
---------------------------------

How is the pipeline actually run? Discuss structure: there is a module file for each broad step. Within that module is a function meant to act as the main script to run that step.

installation instructions :ref:`here <target_to_installation>`
Structure:

- main reduction running script that compiles all of the different steps and runs them.
- The individual modules are a combination of constituent functions and a "wrapper" function that combines the constituent functions to run the step.
- Options are governed by a config file described below.
- Detail the directory structure

There is an example reduction script include in the package called "reduce_tull_spectra.py".

If running reduction with default settings on default setup, can easily run!

Navigate to the directory containing the coude data. For example, a directory called "coude_data" with subdirectories containing the raw data for individual nights following a "YYYYMMDD" naming convention.

Copy in the config file example provided with the package and edit the ``working_dir`` variable to by the date you want to reduce (and change any other options you may want to).

Then can call the base package reduction run script from the command line:

python -m tull_coude_reduction.reduce_tull_spectra config_file_name.yml

And the pipeline will run!

.. _target_to_config_description:

The reduction configuration file
--------------------------------

Below are tables defining the options in each of the main *config* file sections, which govern how the reduction modules are executed.

Each table has three columns: the *config* keyword, a description, and the default value for our base use case with the Tull coudé spectrograph for the echelle setup used by Professor Adam Kraus's group at UT Austin.

:blue:`Every module's config block begins with a` ``do_step`` :blue:`flag that determines if a module is executed.`

General reduction run settings
++++++++++++++++++++++++++++++

*Config* entries for the top level script to run the reduction and analysis pipeline. Mostly concerning path definition, and some general flags for performing certain actions.

================================= ============================================================================== ===================================================================================
**Config Key**  				  **Description**															     **Default value**
--------------------------------- ------------------------------------------------------------------------------ -----------------------------------------------------------------------------------
``do_all_steps``    			  Flag that, if True, overrides all individual modules ``do_step`` flags         False
``remove_images_at_end``    	  Flag to remove processed image files after pipeline completion                 False
``working_dir``                   Path to the night's data directory (typically folder with date name" YYYYMMDD) '/path/to/data/YYYYMMDD/'
``reduction_dir``                 Name of the sub-directory in ``working_dir`` to hold reduction output          'reduction/'
``sub_dir_list``                  List of sub-directories to create in the night's reduction directory           [ 'cals', 'trace', 'object_files', 'spectrum_files', 'wavecal', 'radial_velocity' ]
``header_info_file_name``         File name for the header information CSV file                                  'header_info.csv'
================================= ============================================================================== ===================================================================================

CCD calibration files module
++++++++++++++++++++++++++++

*Config* entries for the :py:meth:`ccd_calibrations <modules.ccd_calibrations>` module.

================================= =========================================================================== ==============================
``bias_image_type``    			  Header image type keyword value for the bias frames						  'zero'
``flat_image_type``    			  Header image type keyword value for the flat frames						  'flat'
``flat_object_names_to_discard``  Header object keyword values to exclude from the creation of the flat field [ 'FF integration time test' ]
``bias_bpm_percentile``    		  Bias flux percentile level above which to mark a bad pixel                  99.9
``flat_field_bpm_limit``          Flat field response value below which to mark a bad pixel                   1.0e-4
================================= =========================================================================== ==============================

Image processing module
+++++++++++++++++++++++

*Config* entries for the :py:meth:`image_processing <modules.image_processing>` module.

================================= ==================================================================== ==============================
``arc_lamp_image_type``    	      Header image type keyword value for arc lamp spectra                 'comp'
``valid_arc_lamp_object_names``   Header object keyword values to include for arc lamp frames (a list) [ 'thar', 'a' ]
``invalid_science_object_names``  Header object keyword values to exclude as science frames (a list)   [ 'test' ]
``cosmic_subtract``               Flag for whether or not to perform cosmic ray subtraction            True
``cosmic_subtract_niter``         Number of iterations for cosmic ray subtraction                      4
================================= ==================================================================== ==============================

Order tracing module
++++++++++++++++++++

*Config* entries for the :py:meth:`trace_echelle <modules.trace_echelle>` module.

============================== ======================================================================================================================= ==============================
``order_xdisp_trace_width``    Cross dispersion pixel height of the slit                                                                               20
``order_start_index``          The index at which to find the starting location of the order traces                                                    -33
``order_center_method``        Algorithm to use for finding the initial locations of the order traces. Must be 'peak_find' or 'gradient_threshold'     'peak_find'
``trace_poly_degree``          Polynomial degree to fit to the trace values                                                                            2
``trace_poly_fit_start_index`` The starting dispersion pixel to use for fitting the trace polynomial                                                   512
``number_of_orders``           Number of orders to trace. If fewer are found, the trace is extended                                                    58
============================== ======================================================================================================================= ==============================

Spectral extraction module
++++++++++++++++++++++++++

*Config* entries for the :py:meth:`extract_spectrum <modules.extract_spectrum>` module.

=============================== ======================================================================================================================= ==============================
``reverse_traced_orders``       Flag for whether or not the order direction needs to be reversed (to match increasing wavelength order)                 True
``order_xdisp_width_extract``   Cross dispersion pixel width of an order for extraction                                                                 -33
``science_extract_type``        Extraction algorithm to use for on-sky science observations. Must be 'optimal_extractin' or 'sum_extraction'            'optimal_extraction'
``lamp_extract_type``           Extraction algorithm to use for arc lamp observations. Must be 'optimal_extractin' or 'sum_extraction'                  'sum_extraction'
``science_background_subtract`` Background subtraction method to use for on-sky science observations. Must be 'fixed' or 'fit'                          'fixed'
``lamp_background_subtract``    Background subtraction method to use for arc lamp observations. Must be 'subtract' or 'none'                            'none'
=============================== ======================================================================================================================= ==============================

Wavelength calibration module
+++++++++++++++++++++++++++++

*Config* entries for the :py:meth:`wavelength_solve_and_calibrate <modules.wavelength_solve_and_calibrate>` module.

================================== ======================================================================================================================= ==============================
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

*Config* entries for the :py:meth:`continuum_fit <modules.continuum_fit>` module.

=============================== ======================================================================================================================= ==============================
``lower_sigma_reject``          Number of standard deviations to reject flux values below the continuum fit                                             2
``upper_sigma_reject``          Number of standard deviations to reject flux values above the continuum fit                                             5
``num_spectrum_chunks``         Number of chunks to break the spectrum into -- sets the number of spline knots                                          6
``min_knot_spacing``            Minimum knot spacing in wavelength (angstrom)                                                                           15
=============================== ======================================================================================================================= ==============================

Radial velocity module
++++++++++++++++++++++

*Config* entries for the :py:meth:`radial_velocity <modules.radial_velocity>` module.

=============================== ======================================================================================================================= ==============================
``orders_to_use_file_name``     File name of the pre-defined orders to use for the broadening function computation                                      'good_orders_tull_coude.csv'
``template_file_name``          File name of the spectrum template to use for broadening function computation                                           'phoenix_t5500_g4.5_m0.0.csv'
``bf_velocity_span``            Velocity span to compute the broadening function over (in km/s)                                                         300
``bf_smooth_res``               Resolution to smooth the broadening function to for fitting (normally instrument resolution)                            60000
``n_bootstrap_samples``         Number of bootstrap samples for brodening function comination to measure the radial velocity                            2500
=============================== ======================================================================================================================= ==============================

.. _target_to_installation:

Installation
------------

The pipeline has not yet been formally released, published on PyPI, and is not pip installable.

Regardless, we  would recommend that any users install a development version of the pipeline because it is still under active development. While the pipeline does currently run in full, it is not yet thoroughly tested enough to warrant a stable v1.0 release. There is also some functionality yet to be added (see `this project <https://github.com/users/dkrolikowski/projects/1>`_ in the pipeline's GitHub repository).

To install the pipeline follow these steps:

We recommend creating a separate python environment to contain the pipeline dependencies that require specific versions. :blue:`Note that the pipeline was developed using Python 3.9`. You can create this with either a python virtual environment or conda. While we name the environment ``pyenv_tull_reduce`` here, feel free to replace it. If using a python virtual environment, first navigate to the place you want to put the environment directory

Follow this code if using a virtual environment, after first navigating to the location you want to place the environment directory: ::

	python -m venv pyenv_tull_reduce
	source pyenv_tull_reduce/bin/activate
	pip install --upgrade setuptools wheel pip

Follow this for using conda: ::

	conda env create -n pyenv_tull_reduce python=3.9
	conda activate pyenv_tull_reduce
	conda update pip setuptools wheel

Note that we also update packages used for installation after we activate the environment.

Then, navigate to the directory you would like to contain the reduction pipeline repository (for example ``~/codes/``). Here we will clone the pipeline GitHub repository and install it. The repository will contain the pipeline reduction running script, the modules, and the needed reference data files. Follow this code to clone the repository and install using pip with the editable flag to install the code in place (so that it can be modified as development continues): ::

	git clone git@github.com:dkrolikowski/tull_coude_reduction.git
	pip install -e tull_coude_reduction

With that, you should have access to the reduction modules and can invoke the reduction script included in the ``tull_coude_reduction`` package.

.. note::

	There are a handful of excess dependencies included in the setup.py file to ensure that there are no missing recursive dependencies.

	In testing, issues like this mostly cropped up with ``saphires``. For example, ``numpy`` is set to version 1.23.5 because ``saphires`` still usings numpy type alises (e.g. ``np.float``). 

	The ``urllib3`` version is set to 1.26.15 due to issues when testing installation of ``barycorrpy``.

.. warning::

	We encountered a failure when testing an installation on the University of Arizona HPC. ``saphires`` use of PyQt5 was not valid, perhaps from an issue with using a ``qt`` backend. This will be investigate in the future and we will coordinate with the ``saphires`` package to fix any bugs there related to this. For now, just keep it in mind in case you encounter issues installing.



