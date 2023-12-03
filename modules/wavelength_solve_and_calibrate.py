""" Functions to wavelength calibrate the spectra. Finds lines in arc lamp spectra (ThAr), matches them to a line list based on a preliminary input wavelength solution, and then iteratively fits a polynomial solution.

Created by DMK on 11/11/2023

Last updated by DMK on 11/13/2023
"""

##### Imports
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from astropy.io import fits
from datetime import datetime
from matplotlib.backends.backend_pdf import PdfPages
from scipy import optimize, signal, stats, interpolate

import os

from tull_coude_reduction.modules import reduction_utils

##### Functions

def order_offset_with_wave_sol_guess( prelim_wavelengths, obs_flux, ref_wavelength, ref_flux, order_use_range = [ 12, 40 ], offset_test_radius = 10 ):
    """ Function to determine if there is an order offset between the wavelength solution guess and the observed spectrum.
    It uses a reference arc lamp spectrum to calculate residuals with the observed lamp spectrum while shifting the wavelength solution guess to determine the best fit order offset.
    It primarily works if the wavelength solution guess is already pretty good, and is used mostly for the wavelength solution guess having more orders than extracted (padding the calibration).

    Parameters
    ----------
    prelim_wavelengths : array
        The preliminary wavelength solution guess, with shape (number of orders, number of pixels). There can be a different number of orders than the extracted spectrum.
    obs_flux : array
        The observed arc lamp flux spectrum, with shape (number of extracted orders, number of pixel).
    ref_wavelength : array
        The wavelength array for the arc lamp reference spectrum.
    ref_flux : array
        The flux array for the arc lamp reference spectrum.
    order_use_range : list of int, optional
        List of minimum and maximum orders to calculate the order offset for. The default is [ 12, 40 ].
    offset_test_radius : int, optional
        The order offset +/- each order to calculate residuals for. The default is 10.

    Returns
    -------
    order_offset : int
        The minimum residual order offset between the wavelength solution guess and observed spectrum.
    """
    
    # The array of orders to calculate the best fit offset for
    orders_to_use = np.arange( *order_use_range )
    
    # The offsets to test
    order_offset_test_arr = np.arange( -offset_test_radius, offset_test_radius + 1 )
    
    # Offset to hold the offset with the minimum residual across test orders
    offset_with_min_resid = np.full( orders_to_use.size, np.nan )
    
    # Go through each of the orders to test
    for i_order, order in enumerate( orders_to_use ):
    
        sum_abs_residuals = np.full( order_offset_test_arr.size, np.nan )
    
        # Go through each offset
        for i_offset, offset in enumerate( order_offset_test_arr ):
    
            # "shift" the wavelength solution guess by the order offset amount
            test_wave = prelim_wavelengths[order+offset]
    
            # Interpolate the reference arc spectrum onto the wavelength solution guess
            ref_flux_interp = np.interp( test_wave, ref_wavelength, ref_flux )
    
            # Calculate the residuals between the reference and data spectrum -- normalize each by the median
            residuals = obs_flux[order] / np.nanmedian( obs_flux[order] ) - ref_flux_interp / np.nanmedian( ref_flux_interp )
                
            # Calculate and output the summed absolute residuals
            sum_abs_residuals[i_offset] = np.nansum( np.abs( residuals ) )
    
        # Find the offset with the minimum absolute residuals for this order
        offset_with_min_resid[i_order] = order_offset_test_arr[np.nanargmin( sum_abs_residuals )]

    # Find the minimum residual offset that is most common across all test orders
    unique_order_offsets, order_offset_counts = np.unique( offset_with_min_resid, return_counts = True )
    order_offset = int( unique_order_offsets[np.argmax(order_offset_counts)] )
    
    return order_offset

def get_flux_mad_from_spectral_chunks( obs_flux, chunk_size = 50 ):
    """ Function to get a noise estimate for the arc lamp spectrum by looking at the MAD of chunks of the spectrum.
    Designed to reduce the influence of bleeding over of saturated Ar lines (and oxide bands to a lesser extent)

    Parameters
    ----------
    flux : array
        The observed spectrum flux array.
    chunk_size : int, optional
        The number of pixels in the spectral chunks to calculate the MAD within. The default is 50.

    Returns
    -------
    median_chunk_mad : float
        The median MAD calculated across all the spectral chunks.
    """
    
    # The number of chunks from the input chunk size in pixels
    number_of_chunks = obs_flux.size // chunk_size
    
    # The array to hold the MAD for each of the spectral chunks
    mad_arr = np.full( number_of_chunks, np.nan )
    
    # Go through each of the chunks
    for i in range( number_of_chunks ):
    
        # The indices of the spectrum within this chunk
        ind_low, ind_high = chunk_size * i, chunk_size * ( i + 1 )
    
        # Calculate the scaled MAD
        mad_arr[i] = stats.median_abs_deviation( obs_flux[ind_low:ind_high], nan_policy = 'omit', scale = 'normal' )
    
    # Adopt the median of the chunk MAD array
    median_chunk_mad = np.nanmedian( mad_arr )

    return median_chunk_mad

def find_arc_lamp_line_pixel_centers( flux, config ):
    """ Function to find lines in an input lamp spectrum using scipy.signal's find_peaks algorithm, with options set in the config file.

    Parameters
    ----------
    flux : array
        The flux spectrum array.
    config : dict
        The overall config file defined using YAML with all parameters for running the reduction and analysis pipeline.

    Returns
    -------
    peak_pixels_initial_fit : array
        The array of pixel centroids found as lines in the input flux spectrum.
    """
    
    ### First prepare from input flux spectrum: normalize it and get noise estimate
    
    # Normalize the flux by its median
    flux = flux / np.nanmedian( flux )
    
    # Now get an estimate of the noise for peak finding threshold. Two current options: one is just from the MAD of the flux spectrum, the other divides the spectrum into chunks and takes the median of each chunk's MAD.
    if config['wavecal']['peak_threshold_mad_method'] == 'full_spectrum':
        
        mad_peak_threshold = stats.median_abs_deviation( flux, nan_policy = 'omit', scale = 'normal' )

    elif config['wavecal']['peak_threshold_mad_method'] == 'chunk_spectrum':
        
        mad_peak_threshold = get_flux_mad_from_spectral_chunks( flux )
        
    ### Now find the peaks using scipy signal's find peak algorithm, with the prominence threshold based on the "noise" estimate of the spectrum
    
    peak_pixels_initial, peak_properties = signal.find_peaks( flux, distance = config['wavecal']['lamp_line_min_separation_pix'], width = config['wavecal']['lamp_line_pix_width_limits'], prominence = config['wavecal']['lamp_line_peak_threshold_sigma'] * mad_peak_threshold )
    
    # Now go through each of the lines found and fit with a Gaussian to get the decimal pixel centroid
    
    peak_pixels_initial_fit = np.full( peak_pixels_initial.size, np.nan )
    
    for i_peak, peak in enumerate( peak_pixels_initial ):
        
        # Get the fit location: within +/- the config 'lamp_line_min_separation_pix' of the peak centroid
        fit_range = np.arange( peak - config['wavecal']['lamp_line_min_separation_pix'], peak + config['wavecal']['lamp_line_min_separation_pix'] + 1 )
        # Get rid of any nans!
        fit_range = np.intersect1d( fit_range, np.where( np.isfinite( flux ) )[0] )

        # The initial guess for 1D gaussian line parameters: the flux at the found centroid, the found int centroid, the minimum of the config 'lamp_line_pix_width_limits', and 0
        p_guess = [ flux[peak], peak, min( config['wavecal']['lamp_line_pix_width_limits'] ), 0 ]
        
        # Fit!
        try:
            line_fit, _ = optimize.curve_fit( reduction_utils.gaussian_1d, fit_range, flux[fit_range], p0 = p_guess )
        except:
            continue
        
        # Output the new fit centroid
        peak_pixels_initial_fit[i_peak] = line_fit[1]
        
    # Get rid of any nans!
    peak_pixels_initial_fit = peak_pixels_initial_fit[np.isfinite( peak_pixels_initial_fit )]
        
    return peak_pixels_initial_fit

def fit_wavelength_solution( pixel_centroids, prelim_wavelengths, line_list_wavelengths, config ):
    """ Function to fit a polynomial wavelength solution, iterating based on sigma rejection with the fit.

    Parameters
    ----------
    pixel_centroids : array
        The array of pixel centroid locations (dependent variable for the wavelength solution).
    prelim_wavelengths : array
        A preliminary guess for the wavelength solution (shape same as the output spectrum).
    line_list_wavelengths : array
        The array of line list wavelengths for the arc lamp used to match to.
    config : dict
        The overall config file defined using YAML with all parameters for running the reduction and analysis pipeline.

    Returns
    -------
    wave_poly_fit : array
        The final adopted polynomial fit coefficients for the wavelength solution (wavelength as a function of pixel). Size is the polynomial degree + 1.
    line_centroid_record : dict
        A dictionary containing a record of the fitting iterations. Each entry is a list of arrays, with each entry being a subsequent fitting iteration.
        Contains keys:
        'pixel': the pixel centroids used in the fit iteration.
        'wavelength': the corresponding wavelength centroids use in the fit iteration.
        'vel_resid': the residuals between the input wavelength centroids and polynomial fit values in velocity (km/s)
        'poly_coeffs': that iteration's polynomial fit coefficients. The last entry is the same as the output 'wave_poly_fit'
    """
    
    # Turn the pixel centroids into wavelength centroids based on the initial wavelength solution guess
    wavelength_centroids_init = np.interp( pixel_centroids, np.arange( prelim_wavelengths.size ), prelim_wavelengths )
    
    ### Now get the line list wavelength for the found peaks
    
    # Get the absolute difference between the initial fit wavelength centroids and the line list
    wave_diff_with_line_list = np.abs( wavelength_centroids_init[:,np.newaxis] - line_list_wavelengths )
    
    # Get the line list wavelengths for the closest match between observation and line list
    closest_line_wavelengths = line_list_wavelengths[np.nanargmin( wave_diff_with_line_list, axis =1 )]
        
    # Only use lines that have a match in the line list within some threshold of wavelength
    line_use_mask = np.abs( wavelength_centroids_init - closest_line_wavelengths ) < config['wavecal']['max_wave_diff_with_list']
    
    # Dictionary to hold the sequence of iterative wavelength solution fit results: the input pixel-wavelength pairs, velocity residuals, and polynomial coefficients
    # Initialize the pixel and wavelength lists with the array of pixel-wavelength pairs that pass the line list closesness cut
    line_centroid_record = { 'pixel': [ pixel_centroids[line_use_mask] ], 'wavelength': [ closest_line_wavelengths[line_use_mask] ], 
                            'vel_resid': [], 'poly_coeffs': [] }
    
    ### Now run the loop for iteratively fitting the wavelength solution -- set break if the number of points include is fewer than the polynomial degree + 1
    while line_centroid_record['pixel'][-1].size >= ( config['wavecal']['wave_cal_poly_order'] + 1 ):
        
        # Fit a polynomial to the lines!
        wave_poly_fit = np.polyfit( line_centroid_record['pixel'][-1], line_centroid_record['wavelength'][-1], config['wavecal']['wave_cal_poly_order'] )
        line_centroid_record['poly_coeffs'].append( wave_poly_fit )
        
        # Calculate the fit polynomial at the found line peaks
        wave_poly_vals = np.polyval( wave_poly_fit, line_centroid_record['pixel'][-1] )
    
        # Calculate residuals between the catalogue wavelength and the polynomial fit, in velocity
        velocity_residuals = ( wave_poly_vals - line_centroid_record['wavelength'][-1] ) / line_centroid_record['wavelength'][-1] * 3e5
        line_centroid_record['vel_resid'].append( velocity_residuals )
        
        # Get the MAD of the velocity residuals
        velocity_residual_mad = stats.median_abs_deviation( velocity_residuals, scale = 'normal', nan_policy = 'omit' )
    
        # See how many lines are outside the sigma level defined in the config
        mad_reject = np.abs( velocity_residuals - np.nanmedian( velocity_residuals ) ) > config['wavecal']['vel_resid_sigma_reject'] * velocity_residual_mad
    
        # If there are no more points to reject -- break out of the loop!
        if mad_reject.sum() == 0 or (~mad_reject).sum() < ( config['wavecal']['wave_cal_poly_order'] + 1 ):
            break
        else:
            line_centroid_record['pixel'].append( line_centroid_record['pixel'][-1][~mad_reject] )
            line_centroid_record['wavelength'].append( line_centroid_record['wavelength'][-1][~mad_reject] )
    
    # Return the final polynomial fit coefficients separate from the record dictionary for ease of access
    
    return wave_poly_fit, line_centroid_record

def wavelength_solution_post_process( wave_sol_coeffs ):
    
    # The output set of wavelength solution polynomial coefficients -- will have rejected orders replaced with the fit to the coefficient vs. order
    output_wave_sol_coeffs = wave_sol_coeffs.copy()
    
    # Go through each of the polynomial coefficients
    for i_par in range( wave_sol_coeffs.shape[1] ):
    
        # Sigma reject fit to the polynomial coefficient vs. order
        fit, x_use, y_use = reduction_utils.polynomial_fit_sigma_reject( np.arange( wave_sol_coeffs.shape[0] ), wave_sol_coeffs[:,i_par], 4, 5, 3, return_data = True )
        
        # The orders that were rejected in the fits
        reject_orders = np.setdiff1d( np.arange( wave_sol_coeffs.shape[0] ), x_use )
        
        # Replace those coeffs with the fit values
        output_wave_sol_coeffs[reject_orders,i_par] = np.polyval( fit, reject_orders )

    return output_wave_sol_coeffs

def interpolate_wavelength_solution( jd_to_interpolate, jd_reference, wavelength_solution_reference ):
    """ Function to interpolate the wavelength solutions measured with the reference sources (e.g. arc lamps) to the times of science observations.

    Parameters
    ----------
    jd_to_interpolate : float
        The JD of the observation to interpolate to.
    jd_reference : array
        The array of JDs of the reference source observations.
    wavelength_solution_reference : array
        The wavelength solution for the reference source observations.

    Returns
    -------
    wavelength_solution_interpolate : array
        The interpolated wavelength solution at the input observation time.
    wavelength_solution_flag : str
        A string denoting what wavelength solution was adopted. Either 'LINTERP' to denote that it was linearly interpolated, or 'CLOSEST" to denote that it could not be interpolated and the closest reference solution was adopted.
    """
    
    # First check if the observation is bounded by reference spectra -- if not just adopt the closest wavelength solution and return a flag
    if jd_reference.min() > jd_to_interpolate or jd_reference.max() < jd_to_interpolate:
        
        wavelength_solution_interpolate = wavelength_solution_reference[np.argmin( np.abs( jd_to_interpolate - jd_reference ) )]
        
        wavelength_solution_flag = 'CLOSEST'
    
    # If it is bounded by arcs -- linear interpolation
    else:
        interp_fn = interpolate.interp1d( jd_reference, wavelength_solution_reference, axis = 0, kind = 'linear' )
        
        wavelength_solution_interpolate = interp_fn( jd_to_interpolate ).astype( float )
        
        wavelength_solution_flag = 'LINTERP'
        
    return wavelength_solution_interpolate, wavelength_solution_flag

### Plotting functions

def plot_spectra_zoom_windows( obs_wavelength, obs_flux, ref_wavelength, ref_flux, lines_used_wavelength, file_name, window_size = 10, number_of_subplots = 6 ):
    """ Function to plot zoom-ins of the final adopted wavelength calibrated observed spectrum compared to a reference source spectrum. Each page has multiple panels with small cut-out of the spectrum.
    This is to highlight the quality of the wavelength solution, and see how it compares to the reference for individual lines across the spectrum.

    Parameters
    ----------
    obs_wavelength : array
        The final adopted wavelength solution for the observed spectrum.
    obs_flux : array
        The observed flux array.
    ref_wavelength : array
        The reference source wavelength array.
    ref_flux : array
        The reference source flux array.
    lines_used_wavelength : array
        The wavelengths of the reference lines from the line list used to calibrate the wavelength solution.
    file_name : str
        The plot file name to write to.
    window_size : float, optional
        The wavelength width of each individual spectrum window. The default is 10.
    number_of_subplots : int, optional
        The number of spectrum window panels on each page. The default is 6.

    Returns
    -------
    None.
    """
        
    number_of_pages = int( np.ceil( np.ceil( np.ptp( obs_wavelength ) / window_size ) / number_of_subplots ) )
    
    subplot_wave_start = obs_wavelength.min()
    
    with PdfPages( file_name ) as pdf:
        
        for i_page in range( number_of_pages ):
            
            # fig = plt.figure( figsize = ( 12, 8 ), num = 1, clear = True )
            fig = plt.figure( num = 1, clear = True )
            
            fig.set_size_inches( 12, 8 )
            
            i_subplot = 1
            
            while subplot_wave_start <= obs_wavelength.max() and i_subplot <= number_of_subplots:
                
                subplot_wave_end = subplot_wave_start + 12.0
                
                obs_plot_loc  = np.where( ( obs_wavelength >= subplot_wave_start ) & ( obs_wavelength <= subplot_wave_end ) )[0]
                ref_plot_loc  = np.where( ( ref_wavelength >= subplot_wave_start ) & ( ref_wavelength <= subplot_wave_end ) )[0]
                line_plot_loc = np.where( ( lines_used_wavelength >= subplot_wave_start ) & ( lines_used_wavelength <= subplot_wave_end ) )[0]
                
                fig.add_subplot( 230 + i_subplot )
                
                plt.plot( obs_wavelength[obs_plot_loc], ( obs_flux / np.nanmedian( obs_flux ) )[obs_plot_loc], '#323232', lw = 1.25 )
                
                plt.plot( ref_wavelength[ref_plot_loc], ( ref_flux / np.nanmedian( ref_flux ) )[ref_plot_loc] * 5, '#bf3465', lw = 1 )
    
                for line_wavelength in lines_used_wavelength[line_plot_loc]:
                    plt.axvline( x = line_wavelength, c = '#1c6ccc', lw = 0.75, ls = '--' )
                    
                plt.yticks( [], [] )
                plt.gca().set_yscale( 'log' )
                                
                subplot_wave_start += window_size
                i_subplot += 1
                
            pdf.savefig( bbox_inches = 'tight', pad_inches = 0.05 )
            # fig.clear()
            # plt.close( fig )
                            
    return None

def plot_wavelength_fit_iteration_spectra( fit_record, obs_flux, ref_wavelength, ref_flux, file_name ):
    """ Function to plot the wavelength-calibrated observed spectra for each of the fitting iterations compared to a reference spectrum. This is to highlight the quality of the solution.
    It is a multi-page PDF plot where each page is a separate iteration (working backwards -- the first page is the adopted iteration).

    Parameters
    ----------
    fit_record : dict
        A dictionary containing a record of the lines used and fit coefficients for each of the fitting iterations. This is the output from fit_wavelength_solution.
    obs_flux : array
        The observed flux array.
    ref_wavelength : array
        The reference source wavelength array.
    ref_flux : array
        The reference source flux array.
    file_name : str
        The plot file name to write to.

    Returns
    -------
    None.
    """
    
    ### Wrap everything in one multi-page PDF -- one page for each iteration
    with PdfPages( file_name ) as pdf:
        
        ### Go through each iteration, but reverse them -- the last iteration is plotted on the first page, and so on
        for i_iter in reversed( range( len( fit_record['pixel'] ) ) ):
            
            # Make the wavelength solution for this iteration
            wavelength_solution = np.polyval( fit_record['poly_coeffs'][i_iter], np.arange( obs_flux.size ) )
            
            # Make the figure
            fig = plt.figure( num = 1, clear = True )
            fig.set_size_inches( 12, 6 )
            
            # Plot the spectra! Normalize by the median
            plt.plot( wavelength_solution, obs_flux / np.nanmedian( obs_flux ), '#323232', lw = 1.25, label = 'Data Flux' )
            
            # Plot the reference ThAr spectrum. Normalize by the median and then x10 so it is offset from the data flux
            ref_plot_loc = np.where( ( ref_wavelength >= wavelength_solution.min() - 5 ) & ( ref_wavelength <= wavelength_solution.max() + 5 ) )[0]
            plt.plot( ref_wavelength[ref_plot_loc], ref_flux[ref_plot_loc] / np.nanmedian( ref_flux[ref_plot_loc] ) * 10, '#bf3465', lw = 1.0, label = 'Reference Arc Spectrum' )
            
            # Plot the lines
            for line_wavelength in fit_record['wavelength'][i_iter]:
                plt.axvline( x = line_wavelength, c = '#1c6ccc', lw = 0.75, ls = '--' )
            
            ### Labels and such
            plt.xlabel( 'Wavelength (${\\rm\AA}$)' )
            plt.ylabel( 'Flux' )
            plt.title( 'Lines Used: {}'.format( fit_record['wavelength'][i_iter].size ) )

            plt.legend( fontsize = 'small' )
            
            plt.gca().set_yscale( 'log' )
            
            pdf.savefig( bbox_inches = 'tight', pad_inches = 0.05 )
            # plt.gcf().clear()
            # plt.close( plt.gcf() )
            
    plt.cla()

    return None

def plot_wavelength_fit_iteration_residuals( fit_record, vel_resid_sigma_reject, file_name ):
    """ Function to plot the velocity residuals for the iterative wavelength solution fitting. It is a multi-page PDF plot where each page is a separate iteration (working backwards -- the first page is the adopted iteration).

    Parameters
    ----------
    fit_record : dict
        A dictionary containing a record of the lines used and fit coefficients for each of the fitting iterations. This is the output from fit_wavelength_solution.
    file_name : str
        The plot file name to write to.
    vel_resid_sigma_reject : float
        The sigma-level used for rejection, to plot horizontal lines and highlight points that are rejected in each iteration.

    Returns
    -------
    None.
    """
    
    ### Wrap everything in one multi-page PDF -- one page for each iteration
    with PdfPages( file_name ) as pdf:
        
        ### Go through each iteration, but reverse them -- the last iteration is plotted on the first page, and so on
        for i_iter in reversed( range( len( fit_record['pixel'] ) ) ):
            
            # Make the figure
            fig = plt.figure( num = 1, clear = True )
            fig.set_size_inches( 12, 6 )
            
            # Plot the velocity residuals
            plt.plot( fit_record['wavelength'][i_iter], fit_record['vel_resid'][i_iter], 'o', c = '#dfa5e5', mec = '#323232', mew = 0.5 )
            
            # Get the velocity residual MAD
            velocity_residual_mad = stats.median_abs_deviation( fit_record['vel_resid'][i_iter], scale = 'normal', nan_policy = 'omit' )
            
            # Plot x's over the points that are rejected!
            lines_rejected = np.where( np.abs( fit_record['vel_resid'][i_iter] - np.nanmedian( fit_record['vel_resid'][i_iter] ) ) > vel_resid_sigma_reject * velocity_residual_mad )[0]
    
            if lines_rejected.size > 0:
                plt.plot( fit_record['wavelength'][i_iter][lines_rejected], fit_record['vel_resid'][i_iter][lines_rejected], 'x', c = '#bf3465', mew = 1.25, ms = 8 )
            
            # Plot horizontal lines at the +/- sigma rejection level
            for i in [ -1, 1 ]:
                plt.axhline( y = np.nanmedian( fit_record['vel_resid'][i_iter] ) + i * vel_resid_sigma_reject * velocity_residual_mad, ls = '--', lw = 1.25, c = '#323232' )
                            
            ### Labels and such
            plt.xlabel( 'Wavelength (${\\rm\AA}$)' )
            plt.ylabel( 'Fit Velocity Residual (km/s)' )
            
            plt.title( 'Lines Used: {}, Lines Rejected: {}, Residuals Sigma: {:.2f} km/s, Cut-off Sigma: {}'.format( fit_record['wavelength'][i_iter].size, lines_rejected.size, velocity_residual_mad, vel_resid_sigma_reject ) )
            
            pdf.savefig( bbox_inches = 'tight', pad_inches = 0.05 )
            # plt.gcf().clear()
            # plt.close( plt.gcf() )
                        
    return None

##### Main wrapper script for wavelength calibration

def wavelength_solution( file_indices, header_df, config ):
    """ Main function to run for fitting the wavelength solution of a set of reference source spectra (e.g. arc lamps).
    Each reference source observation is read in, lines are found and matched to a line list, and the wavelength solution is iteratively fit with a polynomial. The wavelength solution is appended to the spectrum FITS file.

    Parameters
    ----------
    file_indices : array
        The file indices (in the header information file) of the reference source observations for wavelength solution fitting.
    header_df : pandas DataFrame
        The compiled information from the file headers.
    config : dict
        The overall config file defined using YAML with all parameters for running the reduction and analysis pipeline.

    Returns
    -------
    None.
    """
    
    ##### First make the wavelength solution!
    
    ### Read things in that we need for fitting the wavelength solution
    
    # The preliminary wavelength solution guess
    wavelength_solution_guess = np.load( os.path.join( config['general']['code_data_dir'], config['wavecal']['wave_sol_guess'] ) )
    
    # The arc lamp line list
    lamp_line_list_file = pd.read_csv( os.path.join( config['general']['code_data_dir'], config['wavecal']['line_list'] ) )
    lamp_line_list = lamp_line_list_file['wavelength'].values
    
    # Read in the wavelength vs flux CSV of the ThAr photron reference spectrum
    arc_ref_spectrum = pd.read_csv( os.path.join( config['general']['code_data_dir'], config['wavecal']['arc_ref_file'] ) )
    
    ### Go through each of the input file indices for the arc lamps to use
    for i_file in file_indices:
        
        file_in = fits.open( os.path.join( config['general']['reduction_dir'], 'spectrum_files', 'tullcoude_{}_spectrum.fits'.format( header_df['file_token'].values[i_file] ) ) )
        
        # Make directory for this frame's wavelength solution
        frame_dir_path = os.path.join( config['general']['reduction_dir'], 'wavecal', header_df['file_token'].values[i_file] )
        
        # Make the sub-directories for each of the types of plots
        for plot_dir_name in [ 'fit_residuals', 'fit_iter_spectra', 'adopted_sol_spectra_zoom' ]:
            os.makedirs( os.path.join( frame_dir_path, plot_dir_name ), exist_ok = True )
                
        ### Get the order offset to be applied to the preliminary solution if config is flagged to
        if config['wavecal']['use_prelim_sol_order_offset']:
            order_offset = order_offset_with_wave_sol_guess( wavelength_solution_guess, file_in['extracted flux'].data, arc_ref_spectrum['wavelength'].values, arc_ref_spectrum['flux'].values )
            
            # Don't calibrate with this file if the offset is bad! If it would be accessing an unallowed index of the wavelength solution guess (negative or larger than shape)
            if order_offset < 0 or ( order_offset + file_in[1].data.shape[0] ) > wavelength_solution_guess.shape[0]:
                continue
        else:
            order_offset = 0
        
        all_orders_wave_sol_poly_coeffs = np.full( ( config['trace']['number_of_orders'], config['wavecal']['wave_cal_poly_order'] + 1 ), np.nan )
        all_orders_wave_sol = np.full( file_in[1].data.shape, np.nan )
        
        # Go order by order and get wavelength solution
        for order in range( config['trace']['number_of_orders'] ):
            
            ### Fit the wavelength solution
            
            # Get the peak centroids
            lamp_line_pixel_centroids = find_arc_lamp_line_pixel_centers( file_in[1].data[order], config )
            
            # Fit the wavelength solution
            wavelength_solution_poly_coeffs, wavelength_solution_fit_record = fit_wavelength_solution( lamp_line_pixel_centroids, wavelength_solution_guess[order+order_offset], lamp_line_list, config )
            
            all_orders_wave_sol_poly_coeffs[order] = wavelength_solution_poly_coeffs
            all_orders_wave_sol[order] = np.polyval( wavelength_solution_poly_coeffs, np.arange( all_orders_wave_sol.shape[1] ) )
                        
            ### Make any plots!

            # The fit velocity residuals
            plot_file_name = os.path.join( frame_dir_path, 'fit_residuals', 'fit_residuals_order_{}.pdf'.format( order ) )
            plot_wavelength_fit_iteration_residuals( wavelength_solution_fit_record, config['wavecal']['vel_resid_sigma_reject'], plot_file_name )
            
            # The spectra for each of the fits
            plot_file_name = os.path.join( frame_dir_path, 'fit_iter_spectra', 'fit_iter_spectra_order_{}.pdf'.format( order ) )
            plot_wavelength_fit_iteration_spectra( wavelength_solution_fit_record, file_in[1].data[order], arc_ref_spectrum['wavelength'].values, arc_ref_spectrum['flux'].values, plot_file_name )
            
            # Plot order spectrum with zoom-in windows to highlight the solution
            plot_file_name = os.path.join( frame_dir_path, 'adopted_sol_spectra_zoom', 'spectra_zoom_order_{}.pdf'.format( order ) )
            plot_spectra_zoom_windows( all_orders_wave_sol[order], file_in[1].data[order], arc_ref_spectrum['wavelength'].values, arc_ref_spectrum['flux'].values, wavelength_solution_fit_record['wavelength'][-1], plot_file_name )
    
        ### Output this frame's wavelength solution!
        
        # Re-build the output file rather than append to the input, in case of multiple runs
        output_file = fits.HDUList( [ file_in[0], file_in['extracted flux'], file_in['extracted flux error'], fits.ImageHDU( all_orders_wave_sol, name = 'wavelength' ) ] )
                
        # Add a header keyword to the primary HDU for the polynomial degree of the wavelength solution
        output_file[0].header['WAVPOLYD'] = ( config['wavecal']['wave_cal_poly_order'], 'Polynomial degree of wavelength solution' )
        
        # Add a history entry to the primary HDU to mark that it is wavelength calibrated
        output_file[0].header['HISTORY'] = 'Spectrum wavelength calibrated on {}'.format( datetime.strftime( datetime.now(), '%Y/%m/%d' ) )

        # Write out the file with the wavelength solution -- overwrite the previous file
        output_file.writeto( os.path.join( config['general']['reduction_dir'], 'spectrum_files', 'tullcoude_{}_spectrum.fits'.format( header_df['file_token'].values[i_file] ) ), overwrite = True )
            
    return None

def wavelength_calibrate( obj_file_indices, ref_file_indices, header_df, config ):
    """ Main function to run to wavelength calibrate the observed science spectra.
    The wavelength solutions fit for each of the reference source observations are interpolated onto the observation times of the other spectra.

    Parameters
    ----------
    obj_file_indices : array
        The file indices (in the header information file) of all observations with extracted spectra (including science frames and wavelength reference sources).
    ref_file_indices : array
        The file indices (in the header information file) of the reference source observations for wavelength solution fitting.
    header_df : pandas DataFrame
        The compiled information from the file headers.
    config : dict
        The overall config file defined using YAML with all parameters for running the reduction and analysis pipeline.

    Returns
    -------
    None.

    """
    ### Read in the arc wavelength solutions to use for interpolation
    
    arc_wavelength_solutions = []
    
    for i_file in ref_file_indices:
        
        file_name = os.path.join( config['general']['reduction_dir'], 'spectrum_files', 'tullcoude_{}_spectrum.fits'.format( header_df['file_token'].values[i_file] ) )
        
        file_in = fits.getdata( file_name, extname = 'wavelength' )
        
        arc_wavelength_solutions.append( file_in )
                
    arc_wavelength_solutions = np.array( arc_wavelength_solutions )
            
    ### Go through each of the input file indices for the arc lamps to use
    for i_file in obj_file_indices:
        
        # Interpolate the wavelength solution
        obj_wavelength_solution, obj_wavelength_solution_flag = interpolate_wavelength_solution( header_df['obs_jd'].values[i_file], header_df['obs_jd'].values[ref_file_indices], arc_wavelength_solutions )
        
        ### Output the wavelength calibrated file
        
        # Read in the spectrum file to append wavelength solution to
        file_name = os.path.join( config['general']['reduction_dir'], 'spectrum_files', 'tullcoude_{}_spectrum.fits'.format( header_df['file_token'].values[i_file] ) )
        file_in   = fits.open( file_name )
                
        # Build the output file
        output_file = fits.HDUList( [ file_in[0], file_in['extracted flux'], file_in['extracted flux error'], fits.ImageHDU( obj_wavelength_solution, name = 'wavelength' ) ] )
        
        # Add a header keyword to the primary HDU for the polynomial degree of the wavelength solution
        output_file[0].header['WAVPOLYD'] = ( config['wavecal']['wave_cal_poly_order'], 'Polynomial degree of wavelength solution' )
        
        # Add a header keyword for the wavelength solution interpolation: whether it is linearly interpolated or the closest reference wavelength solution was adopted
        output_file[0].header['WAVTYPE']  = ( obj_wavelength_solution_flag, 'Wavelength solution type, interpolation or closest adoption')
        
        # Add a history entry to the primary HDU to mark that it is wavelength calibrated
        output_file[0].header['HISTORY'] = 'Spectrum wavelength calibrated on {}'.format( datetime.strftime( datetime.now(), '%Y/%m/%d' ) )

        # Write out the file with the wavelength solution -- overwrite the previous file
        output_file.writeto( os.path.join( config['general']['reduction_dir'], 'spectrum_files', 'tullcoude_{}_spectrum.fits'.format( header_df['file_token'].values[i_file] ) ), overwrite = True )

    return None




