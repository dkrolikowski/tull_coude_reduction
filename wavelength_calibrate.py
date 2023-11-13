""" Functions to wavelength calibrate the spectra. Finds lines in arc lamp spectra (ThAr), matches them to a line list based on a preliminary input wavelength solution, and then iteratively fits a polynomial solution.

Created by DMK on 11/11/2023
Last updated by DMK on 11/11/2023
"""

##### Imports
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from astropy.io import fits
from datetime import datetime
from matplotlib.backends.backend_pdf import PdfPages
from scipy import optimize, signal, stats

import os
import tqdm

import tull_coude_utils

##### Functions

def order_offset_with_wave_sol_guess( prelim_wavelengths, flux, arc_ref_wavelength, arc_ref_flux, order_use_range = [ 12, 40 ], offset_test_radius = 10 ):
    """ Function to determine if there is an order offset between the wavelength solution guess and the observed spectrum.
    It uses a reference arc lamp spectrum to calculate residuals with the observed lamp spectrum while shifting the wavelength solution guess to determine the best fit order offset.
    It primarily works if the wavelength solution guess is already pretty good, and is used mostly for the wavelength solution guess having more orders than extracted (padding the calibration).

    Parameters
    ----------
    prelim_wavelengths : array
        The preliminary wavelength solution guess, with shape (number of orders, number of pixels). There can be a different number of orders than the extracted spectrum.
    flux : array
        The observed arc lamp flux spectrum, with shape (number of extracted orders, number of pixel).
    arc_ref_wavelength : array
        The wavelength array for the arc lamp reference spectrum.
    arc_ref_flux : array
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
            arc_ref_flux_interp = np.interp( test_wave, arc_ref_wavelength, arc_ref_flux )
    
            # Calculate the residuals between the reference and data spectrum -- normalize each by the median
            residuals = flux[order] / np.nanmedian( flux[order] ) - arc_ref_flux_interp / np.nanmedian( arc_ref_flux_interp )
    
            # Calculate and output the summed absolute residuals
            sum_abs_residuals[i_offset] = np.nansum( np.abs( residuals ) )
    
        # Find the offset with the minimum absolute residuals for this order
        offset_with_min_resid[i_order] = order_offset_test_arr[np.nanargmin( sum_abs_residuals )]

    # Find the minimum residual offset that is most common across all test orders
    unique_order_offsets, order_offset_counts = np.unique( offset_with_min_resid, return_counts = True )
    order_offset = int( unique_order_offsets[np.argmax(order_offset_counts)] )
    
    return order_offset

def get_flux_mad_from_spectral_chunks( flux, chunk_size = 50 ):
        
    # The number of chunks from the input chunk size in pixels
    number_of_chunks = flux.size // chunk_size
    
    # The array to hold the MAD for each of the spectral chunks
    mad_arr = np.full( number_of_chunks, np.nan )
    
    # Go through each of the chunks
    for i in range( number_of_chunks ):
    
        # The indices of the spectrum within this chunk
        ind_low, ind_high = chunk_size * i, chunk_size * ( i + 1 )
    
        # Calculate the scaled MAD
        mad_arr[i] = stats.median_abs_deviation( flux[ind_low:ind_high], nan_policy = 'omit', scale = 'normal' )
    
    # Adopt the median of the chunk MAD array
    median_chunk_mad = np.nanmedian( mad_arr )

    return median_chunk_mad

def find_arc_lamp_line_pixel_centers( flux, config ):
    """ Function to find lines in an input arc lamp spectrum using scipy.signal's find_peaks algorithm, with options set in the config file.

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
            line_fit, _ = optimize.curve_fit( tull_coude_utils.gaussian_1d, fit_range, flux[fit_range], p0 = p_guess )
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
        if mad_reject.sum() == 0:
            break
        else:
            line_centroid_record['pixel'].append( line_centroid_record['pixel'][-1][~mad_reject] )
            line_centroid_record['wavelength'].append( line_centroid_record['wavelength'][-1][~mad_reject] )
    
    # Return the final polynomial fit coefficients separate from the record dictionary for ease of access
    
    return wave_poly_fit, line_centroid_record

### Plotting functions

def plot_wavelength_fit_iteration_spectra( fit_record, flux, arc_ref_wavelength, arc_ref_flux, file_name ):
    
    ### Wrap everything in one multi-page PDF -- one page for each iteration
    with PdfPages( file_name ) as pdf:
        
        ### Go through each iteration, but reverse them -- the last iteration is plotted on the first page, and so on
        for i_iter in reversed( range( len( fit_record['pixel'] ) ) ):
            
            # Make the wavelength solution for this iteration
            wavelength_solution = np.polyval( fit_record['poly_coeffs'][i_iter], np.arange( flux.size ) )
            
            # Make the figure
            plt.figure( figsize = ( 12, 6 ) )
            
            # Plot the spectra! Normalize by the median
            plt.plot( wavelength_solution, flux / np.nanmedian( flux ), '#323232', lw = 1.25, label = 'Data Flux' )
            
            # Plot the reference ThAr spectrum. Normalize by the median and then x10 so it is offset from the data flux
            ref_plot_loc = np.where( ( arc_ref_wavelength >= wavelength_solution.min() - 5 ) & ( arc_ref_wavelength <= wavelength_solution.max() + 5 ) )[0]
            plt.plot( arc_ref_wavelength[ref_plot_loc], arc_ref_flux[ref_plot_loc] / np.nanmedian( arc_ref_flux[ref_plot_loc] ) * 10, '#bf3465', lw = 1.0, label = 'Reference Arc Spectrum' )
            
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
            plt.close()

    return None

def plot_wavelength_fit_iteration_residuals( fit_record, file_name, vel_resid_sigma_reject ):
    
    ### Wrap everything in one multi-page PDF -- one page for each iteration
    with PdfPages( file_name ) as pdf:
        
        ### Go through each iteration, but reverse them -- the last iteration is plotted on the first page, and so on
        for i_iter in reversed( range( len( fit_record['pixel'] ) ) ):
            
            # Make the figure
            plt.figure( figsize = ( 12, 6 ) )
            
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
            plt.close()
            
    return None

##### Main wrapper script for wavelength calibration

def wavelength_solution_and_calibrate( arc_file_indices, header_df, config ):
    
    ##### First make the wavelength solution!
    
    ### Read things in that we need for fitting the wavelength solution
    
    # The preliminary wavelength solution guess
    wavelength_solution_guess = np.load( os.path.join( config['paths']['code_dir'], 'data', config['wavecal']['wave_sol_guess'] ) )[1:]
    
    # The arc lamp line list
    lamp_line_list = np.load( os.path.join( config['paths']['code_dir'], 'data', config['wavecal']['line_list'] ) )
    
    # Read in the wavelength vs flux CSV of the ThAr photron reference spectrum
    arc_ref_spectrum = pd.read_csv( os.path.join( config['paths']['code_dir'], 'data', config['wavecal']['arc_ref_file'] ) )
    
    ### Go through each of the input file indices for the arc lamps to use
    for i_file in tqdm.tqdm( arc_file_indices ):
        
        file_in = fits.open( os.path.join( config['paths']['reduction_dir'], 'spectrum_files', 'tullcoude_{}_spectrum.fits'.format( header_df['file_token'].values[i_file] ) ) )
        
        # Make directory for this frame's wavelength solution
        frame_dir_path = os.path.join( config['paths']['reduction_dir'], 'wavecal', header_df['file_token'].values[i_file] )
        
        # Make the sub-directories for each of the types of plots
        for plot_dir_name in [ 'fit_residuals', 'fit_iter_spectra', 'final_sol_spectra_zoom' ]:
            os.makedirs( os.path.join( frame_dir_path, plot_dir_name ), exist_ok = True )
        
        print( 'Wavelength solving frame {}'.format( header_df['file_token'].values[i_file] ) )
        
        all_orders_wave_sol_poly_coeffs = np.full( ( config['trace']['number_of_orders'], config['wavecal']['wave_cal_poly_order'] + 1 ), np.nan )
        all_orders_wave_sol = np.full( file_in[1].data.shape, np.nan )
        
        # Go order by order and get wavelength solution
        for order in range( config['trace']['number_of_orders'] ):
            
            ### Fit the wavelength solution
            
            # Get the peak centroids
            lamp_line_pixel_centroids = find_arc_lamp_line_pixel_centers( file_in[1].data[order], config )
            
            # Fit the wavelength solution
            wavelength_solution_poly_coeffs, wavelength_solution_fit_record = fit_wavelength_solution( lamp_line_pixel_centroids, wavelength_solution_guess[order], lamp_line_list, config )
            
            all_orders_wave_sol_poly_coeffs[order] = wavelength_solution_poly_coeffs
            all_orders_wave_sol[order] = np.polyval( wavelength_solution_poly_coeffs, np.arange( all_orders_wave_sol.shape[1] ) )
                        
            ### Make any plots!

            # The fit velocity residuals
            plot_file_name = os.path.join( frame_dir_path, 'fit_residuals', 'fit_residuals_order_{}.pdf'.format( order ) )
            plot_wavelength_fit_iteration_residuals( wavelength_solution_fit_record, plot_file_name, config['wavecal']['vel_resid_sigma_reject'] )
            
            # The spectra for each of the fits
            plot_file_name = os.path.join( frame_dir_path, 'fit_iter_spectra', 'fit_iter_spectra_order_{}.pdf'.format( order ) )
            plot_wavelength_fit_iteration_spectra( wavelength_solution_fit_record, file_in[1].data[order], arc_ref_spectrum['wavelength'].values, arc_ref_spectrum['flux'].values, plot_file_name )
            
        ### Output this frame's wavelength solution!
        
        # Re-build the output file rather than append to the input, in case of multiple runs
        output_file = fits.HDUList( [ file_in[0], file_in['extracted flux'], file_in['extracted flux error'], fits.ImageHDU( all_orders_wave_sol, name = 'wavelength' ) ] )
                
        # Add a header keyword to the primary HDU for the polynomial degree of the wavelength solution
        output_file[0].header['WAVPOLYD'] = ( config['wavecal']['wave_cal_poly_order'], 'Polynomial degree of wavelength solution' )
        
        # Add a history entry to the primary HDU to mark that it is wavelength calibrated
        output_file[0].header['HISTORY'] = 'Spectrum wavelength calibrated on {}'.format( datetime.strftime( datetime.now(), '%Y/%m/%d' ) )

        # Write out the file with the wavelength solution -- overwrite the previous file
        output_file.writeto( os.path.join( config['paths']['reduction_dir'], 'spectrum_files', 'tullcoude_{}_spectrum.fits'.format( header_df['file_token'].values[i_file] ) ), overwrite = True )
        
    return None






