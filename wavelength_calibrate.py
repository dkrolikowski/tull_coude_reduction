""" Functions to wavelength calibrate the spectra. Finds lines in arc lamp spectra (ThAr), matches them to a line list based on a preliminary input wavelength solution, and then iteratively fits a polynomial solution.

Created by DMK on 11/11/2023
Last updated by DMK on 11/11/2023
"""

##### Imports
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from astropy.io import fits
from scipy import optimize, signal, stats

import glob
import os
import tqdm

import tull_coude_utils

##### Functions

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
    """

    Parameters
    ----------
    pixel_centroids : TYPE
        DESCRIPTION.
    prelim_wavelengths : TYPE
        DESCRIPTION.
    line_list_wavelengths : TYPE
        DESCRIPTION.
    config : TYPE
        DESCRIPTION.

    Returns
    -------
    wave_poly_fit : TYPE
        DESCRIPTION.
    line_centroid_record : TYPE
        DESCRIPTION.

    """
    
    # Turn the pixel centroids into wavelength centroids based on the initial wavelength solution guess
    wavelength_centroids_init = np.interp( pixel_centroids, np.arange( 2048 ), prelim_wavelengths )
    
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

##### Main wrapper script for wavelength calibration

def wavelength_solution_and_calibrate( arc_file_indices, header_df, config ):
    
    ##### First make the wavelength solution!
    
    ### Go through each of the input file indices for the arc lamps to use
    for i_file in arc_file_indices:
        
        file_in = fits.open( os.path.join( config['paths']['reduction_dir'], 'spectrum_files', 'tullcoude_{}_spectrum.fits'.format( header_df['file_token'].values[i_file] ) ) )
        
        # Go order by order and get wavelength solution
        for order in range( config['trace']['number_of_orders'] ):
            
            # Get the peak centroids
            lamp_line_pixel_centroids = find_arc_lamp_line_pixel_centers( file_in[1].data[order], config )
            
            # Fit the wavelength solution
            order_wavelength_solution = fit_wavelength_solution( lamp_line_pixel_centroids, config )
            
            # Make any plots!

    ### Then apply the wavelength solution to all of the extracted spectra
    
    return None